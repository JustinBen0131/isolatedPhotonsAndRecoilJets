#!/usr/bin/env python3
"""Receipt-bound post-capture assembly; no submission, deletion or release.

This finite adapter is shared by local qualification and a future batch worker.
The native output may be worker scratch, but this adapter NEVER deletes it.
An immutable assembly receipt is the commit marker; a ROOT file alone is not
accepted publication. JES/centrality physics acceptance remains separate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
import time

import uproot

from canonical_storage import file_sha256
from qualify_canonical_part import qualify, load_augmentation_config

VERSION = "CanonicalCaptureAssemblyV1"
DATASETS = {
    "pp_data": "isPP", "auau_data": "isAuAu",
    "pp_photonjet_sim": "isSim", "pp_inclusivejet_sim": "isSimInclusive",
    "auau_photonjet_sim": "isSimEmbedded",
    "auau_inclusivejet_sim": "isSimEmbeddedInclusive",
    "pp_di_sim": ("isSim", "isSimInclusive"),
}


def _json(path, locator=None):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("duplicate JSON key: " + key)
            result[key] = value
        return result
    with Path(path).open("rb") as stream:
        if locator is None:
            raw = stream.read(16 * 1024**2 + 1)
        else:
            match = re.fullmatch(r"([0-9]+):([1-9][0-9]*):([0-9a-f]{64})", locator)
            if not match:
                raise ValueError('invalid template slice locator')
            offset, size = map(int, match.groups()[:2])
            if size > 65536 or offset + size > Path(path).stat().st_size:
                raise ValueError('template slice exceeds bounds')
            stream.seek(offset)
            raw = stream.read(size)
            if (hashlib.sha256(raw).hexdigest() != match.group(3)
                    or not raw.endswith(b'\n') or raw.count(b'\n') != 1):
                raise ValueError('template slice hash or framing differs')
    if len(raw) > 16 * 1024**2:
        raise ValueError("receipt/request exceeds finite read bound")
    def finite(token):
        number = float(token)
        if not math.isfinite(number):
            raise ValueError("nonfinite JSON number")
        return number
    value = json.loads(raw, object_pairs_hook=unique, parse_float=finite,
                       parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)))
    if not isinstance(value, dict):
        raise ValueError("receipt/request must be an object")
    return value, hashlib.sha256(raw).hexdigest()


def validate_capture(source, request_path, receipt_path, product):
    """Bind actual source bytes, event domain and successful terminal execution.

    This validates the current V5 qualification producer receipt, not arbitrary
    historical Condor receipts. DI retains its role and exact photon/inclusive
    sample identity; it cannot be relabeled as SI to enter this adapter.
    """
    if product not in DATASETS:
        raise ValueError("product lacks a qualified capture receipt adapter")
    request, request_sha = _json(request_path)
    receipt, receipt_sha = _json(receipt_path)
    if request.get('schema') in ('THE327Schema14SIMProductionManifestV1',
                                'THE327Schema14DATAProductionManifestV1'):
        return validate_batch_capture(source, request, receipt, request_sha, receipt_sha, product)
    canonical = (json.dumps(request, sort_keys=True, separators=(",", ":"),
                            allow_nan=False)+"\n").encode()
    if (request.get("schema") != "CleanroomEventRequestV5"
            or request.get("operation") != "RUN_ONE"
            or receipt.get("schema") != "CleanroomEventReceiptV2"
            or receipt.get("request_sha256") != hashlib.sha256(canonical).hexdigest()):
        raise ValueError("producer request/receipt binding mismatch")
    env = request.get("environment", {})
    datasets = DATASETS[product]
    if isinstance(datasets, str): datasets = (datasets,)
    if env.get("RJ_REPLAY_DATASET") not in datasets:
        raise ValueError("producer dataset contradicts product")
    expected_role = ('DATA' if product.endswith('_data') else
                     'DI' if product=='pp_di_sim' else
                     'EMBEDDED' if product.startswith('auau_') else 'SI')
    if env.get("RJ_REPLAY_SI_DI_ROLE") != expected_role:
        raise ValueError("producer DATA/SI/DI/EMBEDDED role mismatch")
    if expected_role == 'DI':
        sample = ('run28_photonjet20_double' if env['RJ_REPLAY_DATASET']=='isSim'
                  else 'run28_jet12_double')
        if (request.get('release') != 'ana.568' or env.get('RJ_REPLAY_SAMPLE') != sample
                or env.get('RJ_SIM_SAMPLE') != sample
                or env.get('RJ_PPG12_PHOTON_YIELD_DOUBLE') != '1'
                or env.get('RJ_PPG12_PERIOD_STRICT_DI') != '1'
                or set(request.get('input_files',{})) != {'g4hit','truthjet'}):
            raise ValueError('DI capture sample/recipe mismatch')
    if (receipt.get("status") != "PASS_EXECUTED_PENDING_LOCAL_ROOT_VALIDATION"
            or type(receipt.get("returncode")) is not int or receipt["returncode"] != 0
            or receipt.get("no_surviving_owned_processes") is not True):
        raise ValueError("producer has no successful terminal cleanup receipt")
    events = request.get("events")
    if type(events) is not int or events <= 0:
        raise ValueError("positive exact event count required")
    start = request.get("source_manifest", {}).get("input_start")
    if type(start) is not int or start < 0:
        raise ValueError("producer source range is not bound")
    progress = receipt.get("progress", {})
    if (progress.get("phase") != "event_loop_complete" or progress.get("latest_error")
            or any(type(progress.get(k)) is not int or progress[k] != start+events
                   for k in ("processed", "total"))
            or any(type(progress.get(k)) is not int or progress[k] != events
                   for k in ("row_processed", "row_total"))
            or any(type(progress.get(k)) is not int or progress[k] != 0
                   for k in ("remaining", "row_remaining", "deterministic_skip_remaining"))):
        raise ValueError("producer event range did not complete")
    if any(
            type(progress.get(k)) is not int or progress[k] != start
            for k in ("deterministic_skip_processed", "deterministic_skip_total")):
        raise ValueError("producer deterministic skip is not bound")
    library_sha = request.get("library", {}).get("sha256")
    if (not isinstance(library_sha, str) or not re.fullmatch(r"[0-9a-f]{64}", library_sha)
            or receipt.get("library_sha256") != library_sha):
        raise ValueError("producer library binding mismatch")
    source = Path(source)
    source_sha = file_sha256(source)
    if (type(receipt.get("output_bytes")) is not int
            or receipt["output_bytes"] != source.stat().st_size
            or receipt.get("output_sha256") != source_sha):
        raise ValueError("native output bytes differ from producer receipt")
    with uproot.open(source, array_cache=None) as root:
        retained = root["ReplayFoundationV1/RJEventV1"].num_entries
        if expected_role == 'DATA' and env.get('RJ_REQUIRE_ORIGINAL_SOURCE_CURSOR') == '1':
            _validate_data_source_coverage(root, start, events, retained, events-retained)
        elif retained != events:
            raise ValueError("native output event count contradicts completed producer")
    return dict(request_file_sha256=request_sha, receipt_file_sha256=receipt_sha,
                request_identity_sha256=receipt["request_sha256"],
                producer_library_sha256=library_sha, native_sha256=source_sha,
                events=events, source_entry_begin=start, source_entry_end=start+events,
                physics_acceptance=False, input_retirement_authorized=False)


def _validate_data_source_coverage(root, start, events, retained, rejected):
    """The original input range is partitioned exactly into saved and skimmed events."""
    base = 'ReplayFoundationV1/'
    expected = dict(input_events=str(events), upstream_rejected_events=str(rejected),
        source_entry_contract='ORIGINAL_PAIRED_DST_CURSOR_V1',
        upstream_rejection_module='CaloStatusSkimmer', upstream_rejection_return='ABORTEVENT')
    if rejected < 0 or retained+rejected != events:
        raise ValueError('DATA input/rejected count mismatch')
    for key, value in expected.items():
        if root[base+key].member('fTitle') != value:
            raise ValueError('DATA native metadata mismatch: '+key)
    rejected_tree = root[base+'RJUpstreamRejectedEventV1']
    if (rejected_tree.num_entries != rejected
            or root[base+'RJEventV1'].num_entries != retained
            or root[base+'RJJetConstituentV1'].num_entries != 0):
        raise ValueError('DATA retained/rejected/constituent count mismatch')
    seen_entries, seen_physical = set(), set()
    for tree in (root[base+'RJEventV1'], rejected_tree):
        rows = tree.arrays(['source_entry_ordinal','physical_event_sequence'],library='np')
        for entry, physical in zip(rows['source_entry_ordinal'],rows['physical_event_sequence']):
            entry, physical = int(entry), int(physical)
            if (not start <= entry < start+events or physical < 0
                    or entry in seen_entries or physical in seen_physical):
                raise ValueError('DATA retained/rejected identity overlap')
            seen_entries.add(entry);seen_physical.add(physical)
    if len(seen_entries) != events:
        raise ValueError('DATA input coverage mismatch')


def validate_batch_capture(source, manifest, receipt, manifest_sha, receipt_sha, product):
    """Exact DATA/SI/DI/embedded receipt; no fabricated foreground receipt.

    The worker's pre/post runtime checks remain mandatory. This additionally
    checks the terminal progress bytes and actual native provenance before
    admitting conversion. DI requires the same pinned recipe as native capture.
    """
    data = manifest.get('schema') == 'THE327Schema14DATAProductionManifestV1'
    receipt_schema = 'THE327Schema14'+('DATA' if data else 'SIM')+'ProductionRowReceiptV1'
    if (manifest.get('status') != 'PASS_FROZEN_UNSUBMITTED'
            or manifest.get('submission_ready') is not True
            or manifest.get('source_schema_version') != 14
            or receipt.get('schema') != receipt_schema
            or receipt.get('status') != 'MEASURED_TERMINAL'
            or receipt.get('production_manifest_sha256') != manifest_sha
            or type(receipt.get('wrapper_exit_code')) is not int
            or receipt['wrapper_exit_code'] != 0 or receipt.get('replay_complete') != 1):
        raise ValueError('batch producer manifest/terminal binding mismatch')
    sample = manifest.get('samples', {}).get(receipt.get('sample_id'), {})
    system = sample.get('system')
    admitted = ('pp_data', 'auau_data') if data else (
        'pp_photonjet_sim', 'pp_inclusivejet_sim', 'auau_photonjet_sim', 'auau_inclusivejet_sim', 'pp_di_sim')
    datasets = DATASETS[product]
    if isinstance(datasets, str): datasets = (datasets,)
    expected_role = 'DATA' if data else 'DI' if product == 'pp_di_sim' else 'SI' if system == 'pp' else 'EMBEDDED'
    if (product not in admitted
            or sample.get('dataset') not in datasets
            or system != product.split('_')[0]
            or sample.get('si_di_role') != expected_role
            or receipt.get('system') != system
            or receipt.get('period') != sample.get('period')):
        raise ValueError('batch producer sample/role mismatch')
    if expected_role == 'DI':
        runtime = manifest.get('runtime', {})
        env = runtime.get('pp_replay_environment', {})
        if ((sample.get('dataset'), sample.get('sample')) not in {
                ('isSim', 'run28_photonjet20_double'), ('isSimInclusive', 'run28_jet12_double')}
                or runtime.get('release') != 'ana.568' or sample.get('period') != 'RUN28'
                or runtime.get('pp_wrapper_sha256') != '88c06070ed6cf33368c7e6597c44a2ddbb426780a0b7be952ba6b1d70f350c28'
                or env.get('RJ_PPG12_PHOTON_YIELD_DOUBLE') != '1'
                or env.get('RJ_PPG12_PERIOD_STRICT_DI') != '1'
                or env.get('RJ_PPG12_PERIOD') not in {'0mrad', '1p5mrad'}
                or not re.fullmatch('[0-9a-f]{64}', str(env.get('RJ_PPG12_DI_RUNTIME_MANIFEST_SHA256', '')))):
            raise ValueError('batch DI capture recipe mismatch')
    events = receipt.get('planned_event_count')
    start = receipt.get('source_entry_begin')
    retained = receipt.get('actual_event_count')
    rejected = receipt.get('upstream_rejected_events', 0)
    if (type(events) is not int or not 0 < events <= 100000
            or type(start) is not int or start < 0
            or type(receipt.get('expected_event_count')) is not int or receipt['expected_event_count'] != events
            or type(retained) is not int or retained < 0
            or type(receipt.get('processed_events')) is not int or receipt['processed_events'] != retained
            or type(rejected) is not int or rejected < 0 or retained+rejected != events
            or (not data and rejected != 0)):
        raise ValueError('batch producer exact event range mismatch')
    if data and (receipt.get('input_events') != events
            or receipt.get('source_identity_validated') is not True
            or receipt.get('source_entry_contract') != 'ORIGINAL_PAIRED_DST_CURSOR_V1'
            or receipt.get('event_count_contract') != 'EXACT_INPUT_EQUALS_RETAINED_PLUS_DOCUMENTED_UPSTREAM_REJECTED_V2'
            or (system == 'auau' and events > 35000)):
        raise ValueError('batch DATA input/rejection contract mismatch')
    for key in ('tuple_records_sha256', 'tuple_event_plan_sha256'):
        expected = sample.get('tuple_plan_sha256' if key == 'tuple_event_plan_sha256' else key)
        if not isinstance(expected, str) or not re.fullmatch('[0-9a-f]{64}', expected) or receipt.get(key) != expected:
            raise ValueError('batch source plan binding mismatch')
    live = receipt.get('live_progress', {})
    progress, progress_sha = _json(live.get('path', ''))
    for key in ('task_id', 'workstream_id'):
        owner = manifest.get(key)
        if (not isinstance(owner, str) or not re.fullmatch(r'[A-Za-z0-9_.:-]+', owner)
                or progress.get(key) != owner):
            raise ValueError('batch producer progress owner mismatch')
    if (live.get('sha256') != progress_sha or progress.get('row_id') != receipt.get('row_id')
            or progress.get('campaign_tag') != manifest.get('campaign_tag')
            or progress.get('phase') != 'event_loop_complete' or progress.get('latest_error')
            or any(type(progress.get(k)) is not int or progress[k] != start+events
                   for k in ('processed', 'total'))
            or any(type(progress.get(k)) is not int or progress[k] != events
                   for k in ('row_processed', 'row_total'))
            or any(type(progress.get(k)) is not int or progress[k] != 0
                   for k in ('remaining', 'row_remaining', 'deterministic_skip_remaining'))
            or any(type(progress.get(k)) is not int or progress[k] != start
                   for k in ('deterministic_skip_processed', 'deterministic_skip_total'))):
        raise ValueError('batch producer progress mismatch')
    runtime = manifest.get('runtime', {})
    library_sha = runtime.get(system + '_library_sha256')
    if not isinstance(library_sha, str) or not re.fullmatch('[0-9a-f]{64}', library_sha):
        raise ValueError('batch producer library binding absent')
    source = Path(source)
    source_sha = file_sha256(source)
    if (type(receipt.get('scientific_output_size_bytes')) is not int
            or receipt['scientific_output_size_bytes'] != source.stat().st_size
            or receipt.get('scientific_output_sha256') != source_sha):
        raise ValueError('batch native output bytes mismatch')
    with uproot.open(source, array_cache=None) as root:
        base = 'ReplayFoundationV1/'
        expected_meta = dict(rj_replay_schema_version='14', rj_replay_complete='1',
                             processed_events=str(retained), code_sha256=library_sha,
                             config_sha256=runtime.get(system + '_config_sha256'))
        if data:
            expected_meta.update(input_events=str(events), upstream_rejected_events=str(rejected),
                source_entry_contract='ORIGINAL_PAIRED_DST_CURSOR_V1',
                upstream_rejection_module='CaloStatusSkimmer',upstream_rejection_return='ABORTEVENT')
        for key, value in expected_meta.items():
            if not isinstance(value, str) or root[base+key].member('fTitle') != value:
                raise ValueError('batch native metadata mismatch: '+key)
        if (root[base+'RJEventV1'].num_entries != retained
                or root[base+'RJJetConstituentV1'].num_entries != 0):
            raise ValueError('batch native count/constituent mismatch')
        if data:
            _validate_data_source_coverage(root, start, events, retained, rejected)
        sources = root[base+'RJSourceOccurrenceV1'].arrays(
            ['dataset', 'sample', 'period', 'si_di_role', 'source_manifest_sha256'], library='np')
        expected_source = {key:sample[key] for key in ('dataset', 'sample', 'period', 'si_di_role')}
        expected_source['source_manifest_sha256'] = manifest_sha
        if any(len(sources[key]) != 1 or str(sources[key][0]) != value
               for key, value in expected_source.items()):
            raise ValueError('batch native source occurrence mismatch')
    return dict(request_file_sha256=manifest_sha, receipt_file_sha256=receipt_sha,
                request_identity_sha256=manifest_sha, producer_library_sha256=library_sha,
                source_manifest_sha256=manifest_sha,
                native_sha256=source_sha, events=events, source_entry_begin=start,
                source_entry_end=start+events, physics_acceptance=False,
                input_retirement_authorized=False)


def validate_assembled_donor(source, request, producer_receipt, product, jes_config):
    """Reuse an exact saved canonical donor without relabeling it as native.

    This is an assembly provenance check, never final product acceptance. The
    qualifier independently reads the retained content and performs the join.
    """
    path = Path(jes_config['assembly_path']).resolve(strict=True)
    assembly, digest = _json(path)
    native, native_digest = _json(producer_receipt)
    _, request_digest = _json(request)
    binding = assembly.get('producer', {})
    donor_digest = file_sha256(source)
    conversion = assembly.get('conversion', {})
    native_hash_key = ('scientific_output_sha256' if native.get('schema') in (
        'THE327Schema14DATAProductionRowReceiptV1', 'THE327Schema14SIMProductionRowReceiptV1')
        else 'output_sha256')
    if (digest != jes_config['assembly_sha256']
            or assembly.get('schema') != VERSION
            or assembly.get('status') != 'PASS_ASSEMBLY_NOT_RELEASE'
            or assembly.get('product') != product
            or assembly.get('physics_acceptance') is not False
            or assembly.get('source_retirement_authorized') is not False
            or assembly.get('output_sha256') != donor_digest
            or assembly.get('output_bytes') != Path(source).stat().st_size
            or conversion.get('output_sha256') != donor_digest
            or conversion.get('status') != 'PASS_LOCAL_VALUE_READBACK_NOT_RELEASE'
            or binding.get('request_file_sha256') != request_digest
            or binding.get('receipt_file_sha256') != native_digest
            or assembly.get('native_producer_receipt') != native
            or binding.get('native_sha256') != native.get(native_hash_key)
            or binding.get('physics_acceptance') is not False):
        raise ValueError('canonical donor assembly/native provenance mismatch')
    return binding, dict(path=str(path), sha256=digest, output_sha256=donor_digest)


def finalize(source, output, *, request, producer_receipt, receipt_output,
             product, augmentation_config=None, augmentation_template=None, max_seconds=300,
             max_output_bytes=1024**3, max_rows=2_000_000, progress=None,
             retained_augmentation=False, lossless_compression='inherited',
             augmentation_template_locator=None):
    if lossless_compression not in ('inherited', 'zstd5', 'lzma4'):
        raise ValueError('unknown lossless compression policy')
    if (augmentation_config is None) == (augmentation_template is None):
        raise ValueError('exactly one augmentation config or row template is required')
    if augmentation_template_locator is not None and augmentation_template is None:
        raise ValueError('template locator requires a template catalogue')
    if type(retained_augmentation) is not bool or (retained_augmentation and augmentation_template is not None):
        raise ValueError('retained augmentation needs an explicit bound donor config')
    source, output, request, producer_receipt, receipt_output, config_input = (
        Path(p).resolve() for p in
        (source, output, request, producer_receipt, receipt_output,
         augmentation_config if augmentation_config is not None else augmentation_template))
    if len({source, output, request, producer_receipt, receipt_output, config_input}) != 6:
        raise ValueError("assembly input/output paths must be distinct")
    if output.exists() or receipt_output.exists():
        raise FileExistsError("refusing existing assembly output or receipt")
    if (not output.parent.is_dir() or not receipt_output.parent.is_dir()
            or output.parent.stat().st_dev != receipt_output.parent.stat().st_dev):
        raise ValueError("output and commit receipt require existing directories on the same filesystem")
    if (type(max_seconds) not in (int, float) or not math.isfinite(max_seconds)
            or max_seconds <= 0 or type(max_output_bytes) is not int or max_output_bytes <= 0):
        raise ValueError("finite positive assembly bounds required")
    if source.stat().st_size > max_output_bytes:
        raise ValueError("native source exceeds admitted assembly byte bound")
    started = time.monotonic()
    if progress:
        progress(dict(phase="validate_producer", processed=0, total=1, remaining=1))
    # Historical inputs are NOT relabeled as newly produced captures. Validate
    # the actual augmentation donor with its own native request/receipt, while
    # independently binding and retaining the old source's exact bytes.
    producer_source = source
    retained_binding = None
    retained_template = None
    canonical_donor = None
    if augmentation_template is not None:
        proposed, template_sha = _json(config_input, augmentation_template_locator)
        if proposed.get('schema') == 'RetainedAugmentationTemplateV1':
            if (set(proposed) != {'schema','retained_source','augmentation'}
                    or not isinstance(proposed['retained_source'], str)
                    or not Path(proposed['retained_source']).is_absolute()
                    or not isinstance(proposed['augmentation'], dict)):
                raise ValueError('invalid retained augmentation template')
            source = Path(proposed['retained_source']).resolve()
            if source in {producer_source,output,request,producer_receipt,receipt_output,config_input}:
                raise ValueError('retained template source aliases producer/publication/evidence')
            retained_template = proposed['augmentation']
            if (retained_template.get('schema') != 'NativeAugmentationConfigV1'
                    or 'donor' in retained_template
                    or not isinstance(retained_template.get('binding'), dict)
                    or not isinstance(retained_template.get('donor_binding'), dict)
                    or retained_template['donor_binding'].get('retained_file_sha256','absent') is not None
                    or file_sha256(source) != retained_template['binding'].get('retained_file_sha256')):
                raise ValueError('retained template must pin old bytes and defer only new donor identity')
            if source.stat().st_size > max_output_bytes:
                raise ValueError('retained template source exceeds assembly byte bound')
            retained_augmentation = True
    if retained_augmentation and retained_template is None:
        retained_binding, options, _ = load_augmentation_config(config_input)
        if options.get('donor') is None:
            raise ValueError('retained augmentation requires an explicit donor')
        producer_source = options['donor'].resolve()
        if (producer_source in {source, output, request, producer_receipt, receipt_output, config_input}
                or file_sha256(source) != retained_binding.retained_file_sha256):
            raise ValueError('retained source/donor identity mismatch')
        if producer_source.stat().st_size > max_output_bytes:
            raise ValueError('native donor exceeds admitted assembly byte bound')
    if (retained_augmentation and retained_template is None
            and options.get('retained_jes', {}).get('assembly_path')):
        binding, canonical_donor = validate_assembled_donor(
            producer_source, request, producer_receipt, product, options['retained_jes'])
    else:
        binding = validate_capture(producer_source, request, producer_receipt, product)
    if progress:
        progress(dict(phase="producer_validated", processed=1, total=1, remaining=0))
    remaining = max_seconds - (time.monotonic() - started)
    if remaining <= 0:
        raise TimeoutError("capture validation exhausted assembly deadline")
    # Qualifier stages privately and independently verifies values before linking.
    with tempfile.TemporaryDirectory(prefix=".capture-assembly-", dir=output.parent) as tmp:
        stage = Path(tmp) / "part.root"
        config_sha = template_sha if augmentation_template is not None else file_sha256(config_input)
        bound_config = config_input
        if retained_template is not None:
            template = dict(retained_template)
            template['donor'] = str(producer_source)
            template['donor_binding']['retained_file_sha256'] = binding['native_sha256']
            if template['donor_binding'].get('source_manifest_sha256','absent') is None:
                if 'source_manifest_sha256' not in binding:
                    raise ValueError('deferred donor manifest requires validated batch provenance')
                template['donor_binding']['source_manifest_sha256'] = binding['source_manifest_sha256']
            bound_config = Path(tmp)/'augmentation.json'
            bound_config.write_text(json.dumps(template,allow_nan=False)+'\n')
            retained_binding, _, _ = load_augmentation_config(bound_config)
        elif augmentation_template is not None:
            template = proposed
            # Native content is unknown until capture. A batch source-manifest
            # hash is also deferred to avoid a cycle: manifest -> plan ->
            # template -> manifest. Resolve it only from validated batch
            # provenance, not a caller-supplied replacement identity.
            if (template.get('schema') != 'NativeAugmentationConfigV1'
                    or not isinstance(template.get('binding'), dict)
                    or template['binding'].get('retained_file_sha256', 'absent') is not None
                    or 'donor' in template or 'donor_binding' in template):
                raise ValueError('native row template requires a deferred output hash and no donor')
            template['binding']['retained_file_sha256'] = binding['native_sha256']
            if template['binding'].get('source_manifest_sha256', 'absent') is None:
                source_manifest = binding.get('source_manifest_sha256')
                if not isinstance(source_manifest, str) or not re.fullmatch('[0-9a-f]{64}', source_manifest):
                    raise ValueError('deferred source manifest requires validated batch provenance')
                template['binding']['source_manifest_sha256'] = source_manifest
            bound_config = Path(tmp) / 'augmentation.json'
            bound_config.write_text(json.dumps(template, allow_nan=False)+'\n')
        conversion = qualify(source, stage, product=product,
            augmentation_config=bound_config, max_seconds=remaining,
            max_output_bytes=max_output_bytes, max_rows=max_rows, progress=progress,
            recompress_trees=lossless_compression != 'inherited',
            recompression_codec='zstd5' if lossless_compression == 'inherited' else lossless_compression,
            native_reco_capture_domain=1 if product.endswith("_sim") and not retained_augmentation else None,
            native_truth_capture=product.endswith("_sim") and not retained_augmentation,
            capture_binding=None if retained_augmentation else binding,
            **(dict(jet_donor_request=request, jet_donor_receipt=producer_receipt)
               if retained_augmentation and canonical_donor is None else {}))
        producer_digest = canonical_donor['output_sha256'] if canonical_donor else binding['native_sha256']
        if (file_sha256(producer_source) != producer_digest
                or (canonical_donor and file_sha256(canonical_donor['path']) != canonical_donor['sha256'])
                or (retained_binding is not None and file_sha256(source) != retained_binding.retained_file_sha256)
                or file_sha256(request) != binding["request_file_sha256"]
                or file_sha256(producer_receipt) != binding["receipt_file_sha256"]
                or (_json(config_input, augmentation_template_locator)[1]
                    if augmentation_template_locator is not None else file_sha256(config_input)) != config_sha):
            raise ValueError("producer source or evidence changed during assembly")
        if time.monotonic() - started > max_seconds:
            raise TimeoutError("capture assembly deadline exceeded before publication")
        result = dict(schema=VERSION, status="PASS_ASSEMBLY_NOT_RELEASE", product=product,
            producer=binding, conversion=conversion, output_sha256=conversion["output_sha256"],
            output_bytes=conversion["output_bytes"], physics_acceptance=False,
            source_retirement_authorized=False)
        if retained_augmentation:
            from dataclasses import asdict
            result['retained_input'] = asdict(retained_binding)
            result['assembly_mode'] = 'RETAINED_SOURCE_PLUS_NATIVE_AUGMENTATION_AND_JES_DONOR'
            if canonical_donor:
                result['assembly_mode'] = 'RETAINED_SOURCE_PLUS_CANONICAL_AUGMENTATION_AND_JES_DONOR'
                result['canonical_donor'] = canonical_donor
        if augmentation_template is not None:
            result['augmentation_template_sha256'] = config_sha
            if augmentation_template_locator is not None:
                result['augmentation_template_locator'] = augmentation_template_locator
        result['native_producer_receipt'], _ = _json(producer_receipt)
        marker = Path(tmp) / "receipt.json"
        marker.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n")
        # Never replace existing paths. Receipt is last: consumers require BOTH
        # and the matching content hash; a crash-created orphan is not success.
        os.link(stage, output)
        os.link(marker, receipt_output)
    if progress:
        progress(dict(phase="assembly_complete",processed=1,total=1,remaining=0))
    return result


def assembly_progress_identity(request_path, producer_receipt_path, output):
    """Keep assembly ownership separate from the producer's event counters."""
    request, _ = _json(request_path)
    receipt, _ = _json(producer_receipt_path)
    if request.get('schema') == 'CleanroomEventRequestV5':
        owner = request.get('progress_identity', {})
    else:
        owner = request
    identity = {key: owner.get(key) for key in ('task_id', 'workstream_id', 'campaign_tag')}
    identity['row_id'] = receipt.get('row_id') or owner.get('row_id')
    if any(not isinstance(value, str) or not value or len(value) > 256
           or any(c.isspace() for c in value) for value in identity.values()):
        raise ValueError('assembly progress requires exact producer ownership')
    identity['assembly_output'] = str(Path(output).resolve())
    return identity


def assembly_progress_record(value, *, identity, product, source, started, now=None):
    """Count committed source files; detailed table work has separate units.

    A 1000-event input is reused data, never 1000 newly processed DST events.
    Per-table counters can restart without making the file counter regress.
    """
    now = time.monotonic() if now is None else now
    complete = value.get('phase') == 'assembly_complete'
    count = int(complete)
    elapsed = max(0.0, now-started)
    return dict(schema='LongRunningCanaryProgressV1', **identity,
        executor='CANONICAL_CAPTURE_ASSEMBLY',
        phase=value.get('phase', 'value_readback'), product=product,
        source_name=Path(source).name, processed_unit='committed_source_files',
        processed=count, total=1, remaining=1-count,
        row_processed=count, row_total=1, row_remaining=1-count,
        skip_processed=0, skip_total=0, skip_remaining=0,
        deterministic_skip_processed=0, deterministic_skip_total=0,
        deterministic_skip_remaining=0, new_dst_events=0,
        rate=count/elapsed if elapsed else 0.0,
        rate_per_second=count/elapsed if elapsed else 0.0,
        eta_seconds=0.0 if complete else None,
        updated_at=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        latest_error=value.get('latest_error') or '',
        phase_work=dict(value))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("source", "output"):
        p.add_argument(name, type=Path)
    for name in ("request", "producer-receipt", "receipt-output", "progress-output"):
        p.add_argument("--"+name, type=Path, required=True)
    augment = p.add_mutually_exclusive_group(required=True)
    augment.add_argument('--augmentation-config', type=Path)
    augment.add_argument('--augmentation-template', type=Path)
    p.add_argument('--augmentation-template-locator',
                   help='hash-bound byte slice in a sealed JSONL template catalogue')
    p.add_argument("--product", choices=sorted(DATASETS), required=True)
    p.add_argument("--max-seconds", type=float, default=300)
    p.add_argument("--max-output-bytes", type=int, default=1024**3)
    p.add_argument("--max-rows", type=int, default=2_000_000)
    p.add_argument('--lossless-compression', choices=('inherited','zstd5','lzma4'),
                   default='inherited', help='standard ROOT codec; no branch or value changes')
    p.add_argument('--retained-augmentation', action='store_true',
                   help='old retained source plus separately validated new native donor; never rerun JES')
    p.add_argument('--assembly-campaign')
    p.add_argument('--assembly-row')
    a=p.parse_args()
    arguments=vars(a)
    assembly_campaign=arguments.pop('assembly_campaign')
    assembly_row=arguments.pop('assembly_row')
    progress_path=arguments.pop('progress_output').resolve()
    if (progress_path.exists() or not progress_path.parent.is_dir()
            or progress_path in {v.resolve() for v in arguments.values() if isinstance(v,Path)}):
        raise ValueError('progress requires a new, distinct path in an existing directory')
    identity = assembly_progress_identity(a.request, a.producer_receipt, a.output)
    if assembly_campaign is not None or assembly_row is not None:
        if any(not isinstance(v,str) or not re.fullmatch(r'[A-Za-z0-9_.:\-]{1,160}',v)
               for v in (assembly_campaign,assembly_row)):
            raise ValueError('assembly campaign and row must be supplied together')
        identity.update(producer_campaign_tag=identity['campaign_tag'],producer_row_id=identity['row_id'],
                        campaign_tag=assembly_campaign,row_id=assembly_row)
    started = time.monotonic()
    first=True
    last_print=0.0
    def observe(value):
        nonlocal first, last_print
        value=assembly_progress_record(value, identity=identity, product=a.product,
            source=a.source, started=started)
        fd,name=tempfile.mkstemp(prefix='.assembly-progress-',dir=progress_path.parent)
        try:
            with os.fdopen(fd,'w') as stream:
                json.dump(value,stream,allow_nan=False)
            if first:
                os.link(name,progress_path);first=False
            else:
                os.replace(name,progress_path)
        finally:
            Path(name).unlink(missing_ok=True)
        now=time.monotonic()
        if now-last_print>=5 or value['phase'] in ('failed','assembly_complete'):
            print(json.dumps(value),flush=True);last_print=now
    try:
        receipt=finalize(**arguments,progress=observe)
    except Exception as error:
        observe(dict(phase='failed',processed=0,total=1,remaining=1,
                     latest_error=type(error).__name__+': '+str(error)))
        raise
    print(json.dumps(dict(status=receipt["status"], output_bytes=receipt["output_bytes"],
                         output_sha256=receipt["output_sha256"])))


if __name__ == "__main__":
    main()
