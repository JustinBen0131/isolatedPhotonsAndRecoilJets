"""Pure per-source capture templates, bound to an explicit qualified profile.

No DST access or invented file checksums. The caller seals the returned JSON
with the existing production manifest. The finalizer resolves only the newly
produced ROOT hash and cyclic producer-manifest hash from validated receipts.
"""
from dataclasses import asdict
import copy
import hashlib
import json

import augmentation_contract as c


def _definition(name, values):
    """Validate profile definitions using the native reader's typed contract."""
    values = copy.deepcopy(values)
    cls, arrays = {
        'timing': (c.TimingDefinition, ()),
        'mbd': (c.MbdDefinition, ('arm_order',)),
        'npb': (c.NPBDefinition, ('samples', 'feature_names', 'native_bounds')),
        'calo': (c.CaloDefinition, ('tower_nodes',)),
    }[name]
    for key in arrays:
        if key in values:
            if not isinstance(values[key], (list, tuple)):
                raise c.ContractError('definition array required: ' + key)
            values[key] = tuple(c.NativeFloat32Bound(**item) if key == 'native_bounds'
                                else item for item in values[key])
    return cls(**values)


def bind_new_capture_template(profile, *, run, segment, inputs,
                              reconstruction_sha256, calibration_context):
    """Bind a fixed profile to exact DATA/SIM paths and unchanged definitions.

    ``calibration_context`` is an explicit caller-declared payload-basis object,
    not a claim of fresh CDB or physics acceptance. Hash it consistently for
    old/new same-basis joins. Runtime payload verification remains mandatory.
    ``inputs`` preserves the production role order and uses literal NONE for
    absent SIM streams. Current native producer records the identical tuple.
    """
    if (not isinstance(profile,dict) or profile.get('schema')!='NativeAugmentationConfigV1'
            or not isinstance(profile.get('binding'),dict)
            or any(key in profile for key in ('donor','donor_binding','donor_migration_policy','donor_event_superset'))):
        raise c.ContractError('one qualified same-source capture profile required')
    c._sha(reconstruction_sha256,'capture producer')
    if (not isinstance(calibration_context,dict) or not calibration_context
            or calibration_context.get('schema')!='CaptureCalibrationBasisV1'):
        raise c.ContractError('explicit calibration-basis context required')
    raw=json.dumps(calibration_context,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    basis_sha=hashlib.sha256(raw).hexdigest()
    sample=profile['binding'].get('sample')
    roles=('jets','jetcalo') if sample in ('pp_data','auau_data') else ('calo_cluster','g4hits','jets','global','mbd_epd')
    if not isinstance(inputs,dict) or set(inputs)!=set(roles):
        raise c.ContractError('exact ordered production role set required')
    if any(not isinstance(inputs[role],str) for role in roles):
        raise c.ContractError('literal source paths or NONE required')
    tuple_record='\t'.join(inputs[role] for role in roles)+'\n'
    # Internal placeholder only for validating the deferred-output contract;
    # never serialized as an alleged content or producer-manifest checksum.
    binding=c.SourceBinding(run=run,segment=segment,sample=sample,
        inputs=tuple(c.SourceFile(role,inputs[role],None) for role in roles if inputs[role]!='NONE'),
        source_manifest_sha256='0'*64,retained_file_sha256='0'*64,
        reconstruction_sha256=reconstruction_sha256,calibration_sha256=basis_sha,
        input_identity_kind='FROZEN_TUPLE_RECORD_V1',input_tuple_record=tuple_record)
    result=copy.deepcopy(profile)
    result['binding']=asdict(binding)
    result['binding']['retained_file_sha256']=None
    result['binding']['source_manifest_sha256']=None
    semantic_updates = {}
    for name in ('timing','mbd','npb','calo'):
        definition=result.get(name)
        if name=='calo' and definition is None:continue
        if not isinstance(definition,dict):
            raise c.ContractError('required profile definition absent: '+name)
        old = _definition(name, definition)
        definition['reconstruction_sha256']=reconstruction_sha256
        definition['calibration_sha256']=basis_sha
        new = _definition(name, definition)
        semantic_updates[old.semantic_sha256] = new.semantic_sha256
    if 'retained_semantics' in result:
        semantics = result['retained_semantics']
        allowed_fields = set(c.EVENT_FIELDS + c.CALO_FIELDS + c.CANDIDATE_FIELDS)
        if (not isinstance(semantics, dict) or semantics.keys() - allowed_fields
                or any(not isinstance(value, str) or value not in semantic_updates
                       for value in semantics.values())):
            raise c.ContractError('retained semantics must refer to validated profile definitions')
        result['retained_semantics'] = {field: semantic_updates[value]
                                        for field, value in semantics.items()}
    return dict(template=result,calibration_context=copy.deepcopy(calibration_context),
        calibration_context_sha256=basis_sha,input_tuple_sha256=binding.input_tuple_sha256,
        input_identity_kind='FROZEN_TUPLE_RECORD_V1',dst_content_checksums_claimed=False,
        scientific_acceptance=False)
