"""Pure rendering of a sealed, conversion-only capture assembly worker.

No transport, scheduler, DST loop, deletion or caller-provided shell surface.
The worker consumes a saved native ROOT and publishes via the same finalizer
used immediately after new capture. It is deliberately not a producer rerun.
"""
from __future__ import annotations

import re
import shlex
from pathlib import PurePosixPath

MODULES = frozenset(name+'.py' for name in (
    'augmentation_contract', 'canonical_augmentation_io', 'canonical_contract',
    'canonical_repairs_io', 'canonical_storage', 'centrality_replay',
    'finalize_canonical_capture', 'native_augmentation_reader', 'physics_repairs',
    'qa_augmentation', 'qualify_canonical_part', 'canonical_donor_native_view',
    'retained_jes_donor_join', 'jes_payload_adapter'))
SCHEMA = 'RetainedCaptureAssemblyManifestV1'
AUGMENTATION_SCHEMA = 'RetainedAugmentationAssemblyManifestV1'
PRODUCTS = frozenset(('pp_data', 'auau_data', 'pp_photonjet_sim', 'pp_inclusivejet_sim',
                     'auau_photonjet_sim', 'auau_inclusivejet_sim', 'pp_di_sim'))


def remote_path(value):
    if (not isinstance(value, str) or not re.fullmatch(r'/[A-Za-z0-9_./+\-]+', value)
            or str(PurePosixPath(value)) != value or '..' in PurePosixPath(value).parts
            or not value.startswith(('/sphenix/', '/gpfs/mnt/gpfs02/sphenix/'))):
        raise ValueError('exact scientific remote path required')
    return value


def validate(manifest):
    expected = {'schema', 'task_id', 'workstream_id', 'campaign_tag', 'row_id',
                'packet_root', 'product', 'inputs', 'modules', 'output', 'receipt',
                'progress', 'max_seconds', 'max_output_bytes', 'max_rows',
                'new_dst_events', 'qualification_only', 'bulk_submission'}
    if set(manifest) != expected or manifest['schema'] not in (SCHEMA, AUGMENTATION_SCHEMA):
        raise ValueError('closed retained-assembly manifest required')
    for key in ('task_id', 'workstream_id', 'campaign_tag', 'row_id'):
        if not isinstance(manifest[key], str) or not re.fullmatch(r'[A-Za-z0-9_.:\-]{1,160}', manifest[key]):
            raise ValueError('invalid assembly identity '+key)
    if (manifest['new_dst_events'] != 0 or type(manifest['new_dst_events']) is not int
            or manifest['qualification_only'] is not True or manifest['bulk_submission'] is not False
            or manifest['product'] not in PRODUCTS):
        raise ValueError('conversion-only qualification scope required')
    for key, ceiling in [('max_seconds', 300), ('max_output_bytes', 512*1024**2),
                         ('max_rows', 2_000_000)]:
        if type(manifest[key]) is not int or not 1 <= manifest[key] <= ceiling:
            raise ValueError('finite assembly limit differs: '+key)
    packet = remote_path(manifest['packet_root'])
    if manifest['campaign_tag'] not in PurePosixPath(packet).parts:
        raise ValueError('packet root must bind campaign')
    modules = manifest['modules']
    if not isinstance(modules, dict) or set(modules) != MODULES:
        raise ValueError('exact assembler dependency closure required')
    for digest in modules.values():
        if not isinstance(digest, str) or not re.fullmatch('[0-9a-f]{64}', digest):
            raise ValueError('invalid module hash')
    inputs = manifest['inputs']
    input_names = {'source', 'request', 'producer_receipt', 'augmentation'}
    if manifest['schema'] == AUGMENTATION_SCHEMA:
        input_names.add('donor')
        if 'donor_assembly' in inputs or 'jes_payload' in inputs:
            input_names.update(('donor_assembly', 'jes_payload'))
    if not isinstance(inputs, dict) or set(inputs) != input_names:
        raise ValueError('exact source/donor assembly inputs required')
    for name, pin in inputs.items():
        if not isinstance(pin, dict) or set(pin) != {'path', 'sha256', 'size_bytes'}:
            raise ValueError('exact input pin required')
        remote_path(pin['path'])
        maximum = manifest['max_output_bytes'] if name in ('source', 'donor') else 16*1024**2
        if (type(pin['size_bytes']) is not int or not 1 <= pin['size_bytes'] <= maximum
                or not isinstance(pin['sha256'], str) or not re.fullmatch('[0-9a-f]{64}', pin['sha256'])):
            raise ValueError('invalid input size/hash')
    if inputs['augmentation']['path'] != packet+'/augmentation.json':
        raise ValueError('augmentation must be staged with this worker')
    outputs = [remote_path(manifest[key]) for key in ('output', 'receipt', 'progress')]
    paths = [p['path'] for p in inputs.values()]+outputs
    if len(set(paths)) != len(paths) or any(manifest['campaign_tag'] not in PurePosixPath(p).parts for p in outputs):
        raise ValueError('distinct campaign-scoped publication required')
    if any(p.startswith(packet+'/') for p in outputs):
        raise ValueError('publication cannot modify the sealed packet')
    return manifest


def render_worker(manifest):
    m = validate(manifest)
    q = shlex.quote
    packet = m['packet_root']
    lines = ['#!/usr/bin/env bash', 'set -euo pipefail', 'umask 022',
             'cd "${_CONDOR_SCRATCH_DIR:?Condor execute scratch is required}"',
             'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONNOUSERSITE=1']
    for name, digest in sorted(m['modules'].items()):
        path = q(packet+'/'+name)
        lines.extend(['test -f '+path+' && test ! -L '+path,
                      'test "$(sha256sum -- '+path+' | cut -d\' \' -f1)" = '+q(digest)])
    # JSON key ordering is not semantic: sealing and replay must emit identical bytes.
    for name in sorted(m['inputs']):
        pin = m['inputs'][name]
        path = q(pin['path'])
        lines.extend(['test -f '+path+' && test ! -L '+path,
                      'test "$(stat -c %s -- '+path+')" = '+str(pin['size_bytes']),
                      'test "$(sha256sum -- '+path+' | cut -d\' \' -f1)" = '+q(pin['sha256'])])
    for name in ('output', 'receipt', 'progress'):
        lines.append('test ! -e '+q(m[name])+' && test ! -L '+q(m[name]))
    lines += ['set +u', 'unset PYTHONPATH',
              'source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.568', 'set -u']
    args = ['/usr/bin/timeout', '--signal=TERM', '--kill-after=30s', str(m['max_seconds'])+'s',
            'python3', packet+'/finalize_canonical_capture.py', m['inputs']['source']['path'], m['output'],
            '--request', m['inputs']['request']['path'],
            '--producer-receipt', m['inputs']['producer_receipt']['path'],
            '--augmentation-config', m['inputs']['augmentation']['path'],
            '--receipt-output', m['receipt'], '--progress-output', m['progress'],
            '--assembly-campaign', m['campaign_tag'], '--assembly-row', m['row_id'],
            '--product', m['product'], '--max-seconds', str(m['max_seconds']),
            '--max-output-bytes', str(m['max_output_bytes']), '--max-rows', str(m['max_rows'])]
    if m['schema'] == AUGMENTATION_SCHEMA:
        # Config is hash pinned above. The finalizer additionally verifies its
        # donor bytes; bind its selected path to the independently pinned input.
        check = ('import json,sys; c=json.load(open(sys.argv[1])); '
                 'assert c.get("donor")==sys.argv[2], "sealed donor path differs"')
        lines.append(shlex.join(['python3','-c',check,m['inputs']['augmentation']['path'],
                                 m['inputs']['donor']['path']]))
        if 'donor_assembly' in m['inputs']:
            check = ('import json,sys; c=json.load(open(sys.argv[1]))["retained_jes"]; '
                     'assert c["assembly_path"]==sys.argv[2] and c["assembly_sha256"]==sys.argv[3]; '
                     'assert c["payload_path"]==sys.argv[4] and c["payload_sha256"]==sys.argv[5]')
            lines.append(shlex.join(['python3','-c',check,m['inputs']['augmentation']['path'],
                m['inputs']['donor_assembly']['path'],m['inputs']['donor_assembly']['sha256'],
                m['inputs']['jes_payload']['path'],m['inputs']['jes_payload']['sha256']]))
        args.append('--retained-augmentation')
    lines.append('exec '+shlex.join(args))
    return ('\n'.join(lines)+'\n').encode()
