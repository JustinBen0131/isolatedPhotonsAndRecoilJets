"""Read exact lossless canonical augmentation columns through their native names."""
from contextlib import contextmanager
import json
import uproot

NAMES={
    'ReplayFoundationV1/RJSourceOccurrenceV1':'sources',
    'ReplayFoundationV1/RJEventV1':'events',
    'ReplayFoundationV1/RJPhotonCandidateV1':'photons',
}
CONFIG='ReplayFoundationV1/config_sha256'

class NativeView:
    def __init__(self, root, *, include_jets=False):
        self.root=root
        self.names={}
        if 'root_native_storage_manifest' not in root:
            return
        manifest=json.loads(str(root['root_native_storage_manifest']))
        if (manifest.get('contract')!='FixedROOTNativeStorageFoundationV1' or
                manifest.get('table_mapping_version')!='native_analysis_names_v1' or
                manifest.get('physical_narrowing') is not False or
                manifest.get('content_dependent_branch_suppression') is not False):
            raise ValueError('exact lossless native canonical storage contract required')
        names=dict(NAMES)
        if include_jets:
            names['ReplayFoundationV1/RJJetV1']='jets'
        for native,physical in names.items():
            record=manifest.get('trees',{}).get(native,{})
            if (native in root or record.get('output_tree')!=physical or physical not in root or
                    getattr(root[physical],'classname',None)!='TTree' or
                    root[physical].num_entries!=record.get('entries')):
                raise ValueError('canonical donor tree identity/count differs: '+native)
            # The frozen source schema must agree with every stored native
            # column. No renaming, invented field, type conversion or join.
            actual=root[physical].typenames()
            if actual!=record.get('typenames'):
                raise ValueError('canonical donor native branch schema differs: '+native)
            self.names[native]=physical
        record=manifest.get('objects',{}).get(CONFIG,{})
        physical='source_objects/'+CONFIG
        if (CONFIG in root or record.get('output_object')!=physical or physical not in root or
                record.get('classname')!=getattr(root[physical],'classname',None)):
            raise ValueError('canonical donor config provenance differs')
        self.names[CONFIG]=physical

    def __getitem__(self,name):return self.root[self.names.get(name,name)]
    def get(self,name):return self.root.get(self.names.get(name,name))

@contextmanager
def open_native_view(path, *, include_jets=False):
    with uproot.open(path) as root:
        yield NativeView(root, include_jets=include_jets)
