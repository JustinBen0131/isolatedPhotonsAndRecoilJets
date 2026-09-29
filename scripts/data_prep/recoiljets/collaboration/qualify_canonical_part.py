#!/usr/bin/env python3
"""Local-only, immutable native-to-canonical qualification of one bounded part.

This is not a production submitter, accepted release, or output-retirement tool.
It composes the fixed native ROOT store and photon repair tables in one staged
file, independently reads both back, and publishes only to a new local path.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import time

import uproot

import canonical_storage as storage
from canonical_repairs_io import (prepare_photon_tables, write_root_tables, verify_root_tables,
                                  prepare_jet_replacements)
from canonical_contract import PRODUCTS
from canonical_augmentation_io import (prepare_augmentation_tables, write_augmentation_tables,
                                      verify_augmentation_tables)

VERSION = "CanonicalPartLocalQualificationV1"
AUGMENTATION_CONFIG_VERSION = "NativeAugmentationConfigV1"


def tree_storage_inventory(root):
    """Read ROOT metadata only; never scan branch arrays for byte accounting."""
    return {name: {"rows": int(root[name].num_entries),
                   "compressed_bytes": int(root[name].compressed_bytes),
                   "uncompressed_bytes": int(root[name].uncompressed_bytes)}
            for name, kind in root.classnames(recursive=True, cycle=False).items()
            if kind == "TTree"}


@dataclass(frozen=True)
class AugmentationInput:
    report: object
    timing: object
    mbd: object
    npb: object
    capture_receipt: dict | None = None
    calo: object | None = None


def load_augmentation_config(path):
    """Load explicit local bindings, never infer calibration or model authority.

    Definitions use the JSON form of the existing contract dataclasses. The
    captured model_file is a provenance label, not proof of model-file bytes.
    """
    import augmentation_contract as c

    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate augmentation config key: {key}")
            result[key] = value
        return result

    def finite_number(token):
        number = float(token)
        if not math.isfinite(number):
            raise ValueError("nonfinite augmentation config number")
        return number

    path = Path(path).resolve()
    with path.open("rb") as stream:
        payload = stream.read(1024**2 + 1)
    if len(payload) > 1024**2:
        raise ValueError("augmentation config exceeds 1 MiB")
    config = json.loads(payload, object_pairs_hook=pairs,
                        parse_float=finite_number, parse_constant=finite_number)
    required = {"schema", "binding", "timing", "mbd", "npb"}
    allowed = required | {"calo", "donor", "donor_binding", "model_file", "retained_semantics",
                          "donor_migration_policy", "donor_event_superset", "retained_jes"}
    if (not isinstance(config, dict) or not required <= config.keys()
            or config.keys() - allowed or config["schema"] != AUGMENTATION_CONFIG_VERSION):
        raise ValueError("invalid augmentation config schema/fields")

    def construct(cls, value, arrays=(), nested=None):
        if not isinstance(value, dict):
            raise ValueError(f"{cls.__name__} must be an object")
        values = dict(value)
        for key in arrays:
            if key in values:
                if not isinstance(values[key], list):
                    raise ValueError(f"{key} must be a JSON array")
                values[key] = tuple(construct(nested[key], item) if nested and key in nested else item
                                    for item in values[key])
        return cls(**values)  # Existing constructors reject unknown fields and invalid semantics.

    def binding(value):
        return construct(c.SourceBinding, value, ("inputs",), {"inputs": c.SourceFile})

    kwargs = dict(timing=construct(c.TimingDefinition, config["timing"]),
                  mbd=construct(c.MbdDefinition, config["mbd"], ("arm_order",)),
                  npb=construct(c.NPBDefinition, config["npb"],
                      ("samples", "feature_names", "native_bounds"), {"native_bounds": c.NativeFloat32Bound}))
    if config.get("calo") is not None:
        kwargs["calo"] = construct(c.CaloDefinition, config["calo"], ("tower_nodes",))
    if (config.get("donor") is None) != (config.get("donor_binding") is None):
        raise ValueError("donor and donor_binding must be supplied together")
    if config.get("donor") is not None:
        donor = config["donor"]
        if not isinstance(donor, str) or not Path(donor).is_absolute():
            raise ValueError("donor requires an explicit absolute local path")
        kwargs.update(donor=Path(donor), donor_binding=binding(config["donor_binding"]))
    if 'retained_jes' in config:
        if config.get('donor') is None or not isinstance(config['retained_jes'], dict):
            raise ValueError('retained JES requires a donor and explicit configuration')
        kwargs['retained_jes'] = config['retained_jes']
    if 'donor_migration_policy' in config or 'donor_event_superset' in config:
        from native_augmentation_reader import ADDITIVE_CAPTURE_POLICY
        if config.get('donor') is None:
            raise ValueError('donor migration options require an explicit donor')
        if config.get('donor_migration_policy') not in (None, ADDITIVE_CAPTURE_POLICY):
            raise ValueError('unknown donor migration policy')
        if type(config.get('donor_event_superset', False)) is not bool:
            raise ValueError('donor event superset must be boolean')
        kwargs.update(donor_migration_policy=config.get('donor_migration_policy'),
                      donor_event_superset=config.get('donor_event_superset', False))
    if "model_file" in config:
        if not isinstance(config["model_file"], str) or not config["model_file"]:
            raise ValueError("model_file must be a nonempty captured path label")
        kwargs["model_file"] = config["model_file"]
    if "retained_semantics" in config:
        semantics = config["retained_semantics"]
        if not isinstance(semantics, dict):
            raise ValueError("retained_semantics must be an object")
        for key, value in semantics.items():
            if key not in c.EVENT_FIELDS + c.CANDIDATE_FIELDS + c.CALO_FIELDS:
                raise ValueError(f"unknown retained field: {key}")
            c._sha(value, "retained semantic hash")
        kwargs["retained_semantics"] = semantics
    return binding(config["binding"]), kwargs, {
        "path": str(path), "sha256": hashlib.sha256(payload).hexdigest(),
        "schema": AUGMENTATION_CONFIG_VERSION,
        "authority": "CALLER_DECLARED_NOT_RUNTIME_QUALIFICATION"}


def qualify(source, output, *, product, max_seconds=300, max_output_bytes=1024**3,
            max_rows=2_000_000, progress=None, augmentation=None, recompress_trees=False,
            augmentation_config=None, native_reco_capture_domain=None, native_truth_capture=False,
            capture_binding=None, jet_repairs=None, jet_donor_request=None, jet_donor_receipt=None,
            recompression_codec="zstd5"):
    """No-overwrite local composition, with source and repaired values read back.

    ``photonResponseLinks`` defaults to all retained objects, not the signal
    selection. Unknown denominator/capture status deliberately stays unknown.
    Existing jet links are preserved. No JES correction is silently substituted.
    """
    if product not in PRODUCTS:
        raise ValueError(f"unknown product: {product}")
    if type(recompress_trees) is not bool:
        raise ValueError("recompress_trees must be a boolean")
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output or output.exists():
        raise ValueError("qualification requires a new output path")
    if not output.parent.is_dir():
        raise ValueError("explicit existing output directory required")
    if not math.isfinite(max_seconds) or max_seconds <= 0 or type(max_output_bytes) is not int or max_output_bytes <= 0:
        raise ValueError("positive finite runtime and output bounds required")
    if type(max_rows) is not int or max_rows <= 0:
        raise ValueError("positive integer row bound required")
    if capture_binding is not None:
        # Caller verifies producer evidence. This is retained provenance, not
        # acceptance authority supplied through a freely named JSON object.
        required={"request_file_sha256","receipt_file_sha256","request_identity_sha256",
                  "producer_library_sha256","native_sha256","events",
                  "source_entry_begin","source_entry_end","physics_acceptance",
                  "input_retirement_authorized"}
        optional={"source_manifest_sha256"}
        if (not isinstance(capture_binding,dict)
                or not required <= set(capture_binding)
                or set(capture_binding) - required - optional
                or capture_binding['physics_acceptance'] is not False
                or capture_binding['input_retirement_authorized'] is not False):
            raise ValueError('invalid diagnostic capture binding')
        import re
        if any(not isinstance(capture_binding[k],str) or
               not re.fullmatch(r'[0-9a-f]{64}',capture_binding[k])
               for k in capture_binding if k.endswith('_sha256')):
            raise ValueError('invalid capture binding hash')
        if ('source_manifest_sha256' in capture_binding and
                capture_binding['source_manifest_sha256'] != capture_binding['request_identity_sha256']):
            raise ValueError('capture binding source manifest differs from producer identity')
        for k in ('events','source_entry_begin','source_entry_end'):
            if type(capture_binding[k]) is not int or capture_binding[k]<0:
                raise ValueError('invalid capture binding range')
        if (capture_binding['events']<=0 or capture_binding['source_entry_end'] !=
                capture_binding['source_entry_begin']+capture_binding['events']):
            raise ValueError('capture binding range mismatch')
        capture_binding=dict(capture_binding)
    if augmentation is not None and augmentation_config is not None:
        raise ValueError("choose augmentation input or config, not both")
    if (jet_donor_request is None) != (jet_donor_receipt is None):
        raise ValueError('jet donor request and receipt must be supplied together')
    if jet_donor_request is not None and (augmentation_config is None or jet_repairs is not None):
        raise ValueError('native jet donor needs augmentation config and no manual repairs')
    started = time.monotonic()

    def remaining():
        value = max_seconds - (time.monotonic() - started)
        if value <= 0:
            raise TimeoutError("canonical qualification deadline exceeded")
        return value

    simulation = product not in ("pp_data", "auau_data")
    expected_sample = ("auau" if product.startswith("auau_") else "pp") + ("_sim" if simulation else "_data")
    config_receipt = None
    jes_join_receipt = None
    if augmentation_config is not None:
        from native_augmentation_reader import decode_native_augmentation
        binding, options, config_receipt = load_augmentation_config(augmentation_config)
        jes_config = options.pop('retained_jes', None)
        if jes_config is not None and jet_repairs is not None:
            raise ValueError('JES may be supplied by only one repair path')
        if binding.sample != expected_sample:
            raise ValueError("augmentation sample contradicts declared product")
        remaining()
        decoded = decode_native_augmentation(source, binding, max_rows=max_rows, **options)
        augmentation = AugmentationInput(decoded.report, options["timing"], options["mbd"],
            options["npb"], decoded.receipt, options.get("calo"))
        donor_capture = None
        if jet_donor_request is not None:
            from finalize_canonical_capture import validate_capture
            from canonical_repairs_io import native_jet_donor_repairs
            if options.get('donor') is None:
                raise ValueError('JES augmentation requires an explicit distinct native donor')
            donor_capture = validate_capture(options['donor'], jet_donor_request,
                                            jet_donor_receipt, product)
            if jes_config is None:
                jet_repairs = native_jet_donor_repairs(source, options['donor'],
                    source_binding=binding, donor_binding=options['donor_binding'],
                    capture_binding=donor_capture, max_rows=max_rows)
        if jes_config is not None:
            from retained_jes_donor_join import configured_join
            jet_repairs, jes_join_receipt = configured_join(source, source_binding=binding,
                options=options, config=jes_config, capture_binding=donor_capture,
                max_rows=max_rows, max_seconds=remaining(), progress=progress)
        remaining()
    prepared = prepare_photon_tables(source, simulation=simulation, max_rows=max_rows,
                                    native_reco_capture_domain=native_reco_capture_domain,
                                    native_truth_capture=native_truth_capture)
    if augmentation is not None and not isinstance(augmentation, AugmentationInput):
        raise TypeError("augmentation must be an explicit AugmentationInput")
    if augmentation is not None and augmentation.report.source.sample != expected_sample:
        raise ValueError("augmentation sample contradicts declared product")
    augmented = (prepare_augmentation_tables(source, augmentation.report, timing=augmentation.timing,
        mbd=augmentation.mbd, npb=augmentation.npb, max_rows=max_rows,
        capture_receipt=augmentation.capture_receipt, calo=augmentation.calo) if augmentation else None)
    source_hash = prepared.receipt["source_sha256"]
    repaired_jets = (prepare_jet_replacements(source, jet_repairs, max_rows=max_rows)
                     if jet_repairs is not None else None)
    if repaired_jets and repaired_jets.receipt['source_sha256'] != source_hash:
        raise ValueError('JES repair source differs from photon repair source')
    jet_replacements = repaired_jets.replacements if repaired_jets else None
    if capture_binding and capture_binding['native_sha256'] != source_hash:
        raise ValueError('capture binding source hash mismatch')
    local_manifest = {"schema": VERSION, "status": "LOCAL_QUALIFICATION_NOT_RELEASE",
        "product": product, "source_sha256": source_hash,
        "product_authority": "CALLER_DECLARED_REQUIRES_RELEASE_CENSUS",
        "photon_association_tree": "photonAssociations", "photon_response_tree": "photonResponseLinks",
        "response_selection": "all_retained_objects_NOT_PHYSICS_SIGNAL_SELECTION",
        "legacy_photon_links": prepared.receipt["legacy_photon_links"],
        "native_photon_links": prepared.receipt["native_photon_links"],
        "jet_link_policy": "UNCHANGED", "jes_status": "NOT_APPLIED_OR_CERTIFIED_BY_THIS_ADAPTER",
        "augmentation_status": "TIMING_MBD_EVENT_TIMING_NPB_NOT_CERTIFIED",
        "centrality_status": "PRESERVED_NOT_CERTIFIED_BY_THIS_ADAPTER",
        "release_acceptance": False, "source_retirement_authorized": False,
        "implementation_sha256": {p.name: storage.file_sha256(p) for p in (
            Path(__file__), Path(storage.__file__),
            Path(__file__).with_name("canonical_repairs_io.py"),
            Path(__file__).with_name("physics_repairs.py"),
            Path(__file__).with_name("augmentation_contract.py"),
            Path(__file__).with_name("canonical_augmentation_io.py"),
            Path(__file__).with_name("native_augmentation_reader.py"),
            Path(__file__).with_name("canonical_donor_native_view.py"))}}
    local_manifest["augmentation_tables_present"] = augmented is not None
    local_manifest["augmentation_config"] = config_receipt
    if repaired_jets:
        local_manifest['jes_status'] = repaired_jets.receipt['status']
        local_manifest['jes_repair'] = repaired_jets.receipt
    if capture_binding is not None:
        local_manifest['capture_binding']=capture_binding
    extra_trees = ["photonAssociations", "photonResponseLinks"]
    extra_objects = {"photon_repair_manifest": "TObjString", "local_qualification_manifest": "TObjString"}
    if augmented:
        extra_trees.extend(augmented.tables)
        extra_objects["augmentation_manifest"] = "TObjString"
    with tempfile.TemporaryDirectory(prefix=".canonical-qualification-", dir=output.parent) as tmp:
        stage = Path(tmp) / "part.root"
        native_storage = storage.write_root_native(source, stage, require_native=True,
            max_seconds=remaining(), max_output_bytes=max_output_bytes, progress=progress,
            recompress_trees=recompress_trees, drop_jet_mass=True,
            jet_replacements=jet_replacements, recompression_codec=recompression_codec)
        import ROOT
        root = ROOT.TFile.Open(str(stage), "UPDATE")
        if not root or root.IsZombie():
            raise RuntimeError("cannot open private stage")
        try:
            write_root_tables(root, prepared, max_seconds=remaining(), progress=progress)
            if augmented:
                write_augmentation_tables(root, augmented, max_seconds=remaining(), progress=progress)
            root.cd()
            ROOT.TObjString(json.dumps(local_manifest, sort_keys=True)).Write("local_qualification_manifest")
        finally:
            root.Close()
        if stage.stat().st_size > max_output_bytes:
            raise ValueError("final qualification output exceeds byte bound")
        remaining()
        native_readback = storage.validate_root_native(source, stage,
            check_deadline=remaining, progress=progress,
            extra_trees=extra_trees, extra_objects=extra_objects, expected_drop_jet_mass=True,
            expected_jet_replacements=jet_replacements)
        repair_readback = verify_root_tables(stage, prepared)
        augmentation_readback = verify_augmentation_tables(stage, augmented) if augmented else None
        with uproot.open(stage) as readback:
            if json.loads(str(readback["local_qualification_manifest"])) != local_manifest:
                raise ValueError("local qualification manifest readback mismatch")
            table_storage = tree_storage_inventory(readback)
        if storage.file_sha256(source) != source_hash:
            raise ValueError("source changed during qualification")
        if augmentation and augmentation.capture_receipt:
            capture = augmentation.capture_receipt
            if storage.file_sha256(capture["donor_path"]) != capture["donor_sha256"]:
                raise ValueError("native donor changed during qualification")
        if config_receipt and storage.file_sha256(config_receipt["path"]) != config_receipt["sha256"]:
            raise ValueError("augmentation config changed during qualification")
        remaining()
        receipt = {"schema": VERSION, "status": "PASS_LOCAL_VALUE_READBACK_NOT_RELEASE",
            "jes_donor_join": jes_join_receipt,
            "source_sha256": source_hash, "output_sha256": storage.file_sha256(stage),
            "output_bytes": stage.stat().st_size, "product": product,
            "storage": {"measured_at": datetime.now(timezone.utc).isoformat(),
                "source_bytes": source.stat().st_size,
                "native_store_bytes": native_storage["output_bytes"],
                "final_package_bytes": stage.stat().st_size,
                "appended_tables_and_metadata_bytes": stage.stat().st_size - native_storage["output_bytes"],
                "tree_compression": native_storage["tree_compression"],
                "scope": "THIS_EXACT_PART_NOT_CAMPAIGN_FORECAST",
                "augmentation_tables_present": augmented is not None},
            "native_readback": native_readback, "repair_readback": repair_readback,
            "augmentation_readback": augmentation_readback,
            "augmentation_config": config_receipt,
            "photon_repair": prepared.receipt, "science_acceptance": False,
            "source_retirement_authorized": False, "elapsed_seconds": time.monotonic() - started}
        receipt['jes_repair'] = repaired_jets.receipt if repaired_jets else None
        receipt["storage"]["tables"] = table_storage
        if capture_binding is not None:
            receipt['capture_binding']=capture_binding
        receipt["storage"]["tree_compressed_bytes"] = sum(
            row["compressed_bytes"] for row in table_storage.values())
        receipt["storage"]["table_byte_scope"] = (
            "ROOT_TREE_METADATA; file size remains authoritative and also includes non-tree objects and headers")
        # Atomic no-clobber creation on the same filesystem, after every readback.
        remaining()
        os.link(stage, output)
        return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--product", choices=sorted(PRODUCTS), required=True)
    parser.add_argument("--max-seconds", type=float, default=300)
    parser.add_argument("--max-output-bytes", type=int, default=1024**3)
    parser.add_argument("--max-rows", type=int, default=2_000_000)
    parser.add_argument("--native-reco-capture-domain", type=int, choices=[1],
                        help="verify persisted domain-1 reco census; legacy evidence remains unknown")
    parser.add_argument("--augmentation-config", type=Path,
                        help="explicit NativeAugmentationConfigV1 JSON; no donor auto-discovery")
    parser.add_argument('--jet-donor-request', type=Path,
                        help='validated native request for the same augmentation donor; replace JES once')
    parser.add_argument('--jet-donor-receipt', type=Path,
                        help='matching successful native terminal receipt; no old-receipt relabeling')
    parser.add_argument("--native-truth-capture", action="store_true",
                        help="cross-check native truth inventory counts and keyed rows; not signal acceptance")
    parser.add_argument("--recompress-trees", action="store_true",
                        help="opt in to ZSTD5 tree recompression; default retains input compression")
    args = parser.parse_args()
    receipt = qualify(args.source, args.output, product=args.product,
                      max_seconds=args.max_seconds, max_output_bytes=args.max_output_bytes,
                      recompress_trees=args.recompress_trees, max_rows=args.max_rows,
                      augmentation_config=args.augmentation_config,
                      jet_donor_request=args.jet_donor_request, jet_donor_receipt=args.jet_donor_receipt,
                      native_reco_capture_domain=args.native_reco_capture_domain,
                      native_truth_capture=args.native_truth_capture)
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
