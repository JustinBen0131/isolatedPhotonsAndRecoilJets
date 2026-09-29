"""Bounded physical-identity JES join for native and lossless canonical donors.

This composes retained augmentation with saved native JetCalib values. It is
local qualification, not historical calibration acceptance or release authority.
The payload is used as an independent signed-value check; accepted positive pT
is copied from the native donor, never multiplied into an already corrected pT.
"""
import math
import time
import json
from pathlib import Path

import awkward as ak

from canonical_donor_native_view import open_native_view
from canonical_storage import file_sha256
from jes_payload_adapter import float32, historical_route
from physics_repairs import JetPTRepair, RepairError

BASE = 'ReplayFoundationV1/'
KEY_FIELDS = ('algorithm', 'radius', 'input_identity', 'subtraction_identity', 'native_raw_jet_key')


def classify_jet(old, new, signed_pt):
    """Reject changed raw jets and unexplained directions; classify signed JES."""
    for name in (*KEY_FIELDS, 'raw_pt', 'deterministic_order', 'native_jet_key'):
        if old[name] != new[name]:
            raise RepairError('donor changed raw jet: '+name)
    if not math.isfinite(signed_pt) or not math.isfinite(new['corrected_pt']):
        raise RepairError('nonfinite JES result')
    tolerance = max(1e-8, 4 * 2**-23 * abs(new['corrected_pt']))
    if abs(abs(signed_pt) - new['corrected_pt']) > tolerance:
        raise RepairError('signed payload disagrees with native donor pT')
    if any(not math.isfinite(row[name]) for row in (old, new) for name in ('eta', 'phi')):
        raise RepairError('nonfinite jet direction')
    phi_delta = abs(math.remainder(old['phi'] - new['phi'], 2*math.pi))
    if signed_pt < 0:
        if abs(old['eta'] + new['eta']) > 1e-5 or abs(phi_delta - math.pi) > 1e-5:
            raise RepairError('negative native JES lacks the expected vector reversal')
        return 'JES_INVALID_NATIVE_SIGN_PRESERVED'
    if abs(old['eta'] - new['eta']) > 1e-5 or phi_delta > 1e-5:
        raise RepairError('positive JES changed direction')
    return 'CALIBRATED_PT_ONLY'


def inventory(path, max_rows):
    with open_native_view(path, include_jets=True) as root:
        def rows(name):
            tree=root[BASE+name]
            if tree.num_entries > max_rows:
                raise RepairError('retained JES table exceeds row bound: '+name)
            return ak.to_list(tree.arrays(library='ak'))
        occurrences=rows('RJSourceOccurrenceV1')
        if len(occurrences)!=1:
            raise RepairError('one source occurrence required')
        events=rows('RJEventV1'); physical={}; ids={}
        for event in events:
            identity=(event['event_id_hi'],event['event_id_lo'])
            sequence=event['physical_event_sequence']
            if (event['physical_event_sequence_valid'] != 1 or type(sequence) is not int
                    or sequence < 0 or sequence in physical or identity in ids):
                raise RepairError('ambiguous physical event identity')
            ids[identity]=sequence;physical[sequence]=event
        jets={}
        for row in rows('RJJetV1'):
            identity=(row['event_id_hi'],row['event_id_lo'])
            if identity not in ids:
                raise RepairError('jet has no physical event')
            key=(ids[identity], *(row[k] for k in KEY_FIELDS))
            if key in jets:
                raise RepairError('duplicate physical jet identity')
            jets[key]=row
        return occurrences[0],physical,jets


def join(source, donor, *, source_binding, donor_binding, assembly, assembly_sha256,
         evaluator, max_rows=2_000_000, max_seconds=180, progress=None, capture_binding=None):
    if type(max_rows) is not int or max_rows <= 0 or not math.isfinite(max_seconds) or max_seconds <= 0:
        raise RepairError('finite positive JES join bounds required')
    source_sha=file_sha256(source);donor_sha=file_sha256(donor)
    if (source_sha!=source_binding.retained_file_sha256 or donor_sha!=donor_binding.retained_file_sha256
            or source_sha==donor_sha):
        raise RepairError('distinct hash-bound retained source and donor required')
    for name in ('run','segment','sample','inputs','input_identity_kind','calibration_sha256'):
        if getattr(source_binding,name)!=getattr(donor_binding,name):
            raise RepairError('source/donor binding differs: '+name)
    assembly = assembly or {}
    producer=assembly.get('producer',{});native=assembly.get('native_producer_receipt',{})
    if capture_binding is not None:
        if (assembly or assembly_sha256 is not None
                or capture_binding.get('native_sha256') != donor_sha
                or capture_binding.get('producer_library_sha256') != donor_binding.reconstruction_sha256
                or capture_binding.get('physics_acceptance') is not False):
            raise RepairError('exact validated native capture required')
        producer = capture_binding
    elif (assembly.get('schema')!='CanonicalCaptureAssemblyV1'
            or assembly.get('status')!='PASS_ASSEMBLY_NOT_RELEASE'
            or assembly.get('output_sha256')!=donor_sha
            or assembly.get('physics_acceptance') is not False
            or assembly.get('source_retirement_authorized') is not False
            or assembly.get('conversion',{}).get('status')!='PASS_LOCAL_VALUE_READBACK_NOT_RELEASE'
            or producer.get('producer_library_sha256')!=donor_binding.reconstruction_sha256
            or native.get('production_manifest_sha256')!=donor_binding.source_manifest_sha256
            or native.get('wrapper_exit_code')!=0
            or type(native.get('replay_complete')) not in (bool,int)
            or native.get('replay_complete') != 1):
        raise RepairError('exact canonical donor producer/assembly receipt required')
    started=time.monotonic()
    def update(phase, processed, total):
        if progress: progress(dict(phase=phase,processed=processed,total=total,remaining=total-processed))
        if time.monotonic()-started>max_seconds: raise TimeoutError('retained JES join deadline')
    update('index_sources',0,1)
    old_occ,old_events,old=inventory(source,max_rows)
    new_occ,new_events,new=inventory(donor,max_rows)
    for field in ('lane','dataset','sample','run','segment','si_di_role'):
        if old_occ[field]!=new_occ[field]: raise RepairError('source population differs: '+field)
    for occurrence,binding in ((old_occ,source_binding),(new_occ,donor_binding)):
        if (occurrence['run'],occurrence['segment'],occurrence['source_manifest_sha256'])!=(
                binding.run,binding.segment,binding.source_manifest_sha256):
            raise RepairError('stored source occurrence differs from binding')
    if old_events.keys()!=new_events.keys() or old.keys()!=new.keys():
        raise RepairError('retained event or jet population differs')
    for key,event in old_events.items():
        old_z, new_z = event['vertex_z'], new_events[key]['vertex_z']
        if old_z != new_z and not (math.isnan(old_z) and math.isnan(new_z)):
            raise RepairError('physical event vertex differs')
    encoding={k:[old_occ[k],new_occ[k]] for k in ('period','input_uri_hash') if old_occ[k]!=new_occ[k]}
    if 'period' in encoding and encoding['period']!=['0mrad','RUN28']:
        raise RepairError('unqualified period encoding migration')
    repairs=[];invalid=0;functions={};max_error=0.
    total=len(old)
    update('classify_native_jes',0,total)
    for i,(key,row) in enumerate(old.items(),1):
        route=historical_route(row['radius'],row['eta'],old_events[key[0]]['vertex_z'])
        name=route.function_name
        if name not in functions: functions[name]=evaluator._function(name)
        signed=float32(functions[name].Eval(float32(row['raw_pt'])))
        status=classify_jet(row,new[key],signed);bad=status=='JES_INVALID_NATIVE_SIGN_PRESERVED'
        invalid+=bad;max_error=max(max_error,abs(abs(signed)-new[key]['corrected_pt']))
        provenance=dict(source_scope='sha256:'+source_sha,donor_sha256=donor_sha,
            evaluation_semantics='EXACT_SAVED_NATIVE_JETCALIB_OUTPUT',assembly_receipt_sha256=assembly_sha256,
            producer_library_sha256=producer['producer_library_sha256'],payload_sha256=evaluator.payload_sha256,
            payload_function=name,source_occurrence_encoding_differences=encoding,
            historical_non_jes_calibration_equivalence=False,scientific_acceptance=False,
            correction_application='PRESERVE_ORIGINAL_AND_MARK_JES_INVALID' if bad else
            'REPLACE_FROM_NATIVE_DONOR_NEVER_MULTIPLY_CORRECTED_PT')
        if bad: provenance.update(invalid_reason='NEGATIVE_SIGNED_TF1_PT_NATIVE_DIRECTION_REVERSAL',signed_payload_pt=signed)
        repairs.append(JetPTRepair(dict(row) if bad else dict(row,corrected_pt=new[key]['corrected_pt']),
            status,'PRE_JES','UNDEFINED_JES_PRESERVED_ORIGINAL' if bad else 'POST_JES',(),
            'RETAINED_ROWS_ONLY_LOCAL_QUALIFICATION',provenance))
        if i%5000==0: update('classify_native_jes',i,total)
    if file_sha256(source)!=source_sha or file_sha256(donor)!=donor_sha:
        raise RepairError('source changed during JES join')
    update('jes_join_complete',total,total)
    return repairs,dict(status='PASS_EXACT_NATIVE_JES_JOIN_NOT_RELEASE',jets=total,invalid_jes_jets=invalid,
        max_payload_native_pt_error=max_error,events=len(old_events),source_sha256=source_sha,
        donor_sha256=donor_sha,payload_sha256=evaluator.payload_sha256,source_occurrence_encoding_differences=encoding,
        historical_non_jes_calibration_equivalence=False,scientific_acceptance=False)


def configured_join(source, *, source_binding, options, config, capture_binding=None,
                    max_rows=2_000_000, max_seconds=180, progress=None):
    """Production qualifier entry: use explicit hash-bound payload/receipt only."""
    from jes_payload_adapter import LocalRootTF1Evaluator
    required = {'schema', 'payload_path', 'payload_sha256'}
    allowed = required | {'assembly_path', 'assembly_sha256'}
    if (not isinstance(config, dict) or not required <= config.keys()
            or config.keys() - allowed or config['schema'] != 'NativeJESDonorJoinV1'
            or options.get('donor') is None):
        raise RepairError('explicit NativeJESDonorJoinV1 with distinct donor required')
    if config['payload_sha256'] != '76a8788fdb4e0b3859d01361f884e295ea609bb8c1cc188105389214db2de27d':
        raise RepairError('payload differs from selected common pp-v6 nominal')
    payload = Path(config['payload_path'])
    if not payload.is_absolute():
        raise RepairError('absolute pinned JES payload required')
    assembly = None
    assembly_sha = None
    if capture_binding is None:
        if not {'assembly_path','assembly_sha256'} <= config.keys():
            raise RepairError('canonical donor requires exact assembly receipt')
        path = Path(config['assembly_path'])
        if not path.is_absolute() or path.stat().st_size > 1024**2:
            raise RepairError('bounded absolute assembly receipt required')
        assembly_sha = file_sha256(path)
        if assembly_sha != config['assembly_sha256']:
            raise RepairError('assembly receipt hash differs')
        assembly = json.loads(path.read_text())
    elif {'assembly_path','assembly_sha256'} & config.keys():
        raise RepairError('native and canonical donor evidence are mutually exclusive')
    evaluator = LocalRootTF1Evaluator(payload, config['payload_sha256'])
    try:
        result = join(source, options['donor'], source_binding=source_binding,
            donor_binding=options['donor_binding'], assembly=assembly,
            assembly_sha256=assembly_sha, evaluator=evaluator, capture_binding=capture_binding,
            max_rows=max_rows, max_seconds=max_seconds, progress=progress)
        if file_sha256(payload) != config['payload_sha256']:
            raise RepairError('JES payload changed during join')
    finally:
        evaluator.close()
    if assembly is not None and file_sha256(path) != assembly_sha:
        raise RepairError('assembly receipt changed during join')
    return result
