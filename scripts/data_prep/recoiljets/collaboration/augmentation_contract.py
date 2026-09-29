"""Pure keyed donor reconciliation; no ROOT, reconstruction, I/O, or inference.

Callers must verify the bytes behind every hash and extract native physical event
and rebuilt photon keys. A well-formed hash is not proof of those bytes. This
module validates one explicitly bound retained source file at a time; it neither
selects source files nor grants production/scientific acceptance. Every retained
event/candidate gets an overlay row, including absent or incomplete donors.

Audited distinctions (2026-09-20): PhotonClusterBuilder::mean_time averages all
owned calibrated TowerInfo cells, whereas PPG12 plot_cluster_time.C averages only
owned 7x7 non-ZS cells. Raw times are samples; native MbdOut times are ns. No time
difference is computed here: that additionally needs the run-offset convention.
NPB uses its own float/denominator semantics and producer vertex, not necessarily
same-named retained ratios or the event vertex. e22 needs the original owned map;
the ordinary 14-feature tight-ID vector is not the 25-feature NPB vector.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
import re
import struct


NPB_FEATURES = (
    "cluster_Et", "cluster_Eta", "vertexz", "e11_over_e33", "e32_over_e35",
    "e11_over_e22", "e11_over_e13", "e11_over_e15", "e11_over_e17",
    "e11_over_e31", "e11_over_e51", "e11_over_e71", "e22_over_e33",
    "e22_over_e35", "e22_over_e37", "e22_over_e53", "cluster_weta_cogx",
    "cluster_wphi_cogx", "cluster_et1", "cluster_et2", "cluster_et3",
    "cluster_et4", "cluster_w32", "cluster_w52", "cluster_w72",
)
SPLIT_MODEL_SHA256 = "d6086dadac534013cda15cdfb69c1683776d3456d9e439903589653e8ac19eab"
SPLIT_TRAINING_CONFIG_SHA256 = "450039db8afcf6dcb16acc7e738a9ce41d3cd0a476da2f7830aebe539d750c53"
NATIVE_TIMING = "NATIVE_ALL_OWNED_TOWERINFO_E_WEIGHTED_V1"
PPG12_TIMING = "PPG12_OWNED_7X7_NOT_ZS_E_WEIGHTED_V1"
SAMPLES = ("pp_data", "auau_data", "pp_sim", "auau_sim")
EVENT_FIELDS = ("mbd_t0", "mbd_time_south", "mbd_time_north")
CALO_FIELDS = ("emcal_total_energy", "ihcal_total_energy", "ohcal_total_energy", "total_calo_energy")
CANDIDATE_FIELDS = ("cluster_time_raw", "npb_score", "auau_npb_score")


class ContractError(ValueError):
    """A join, provenance claim, or numerical witness is contradictory."""


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ContractError(f"{name} must be nonempty text")


def _sha(value, name):
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ContractError(f"{name} must be a lowercase SHA256")


def _int(value, name, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ContractError(f"invalid integer {name}")


def _finite(value, name):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ContractError(f"{name} must be finite numeric data, not a sentinel state")


def _typed(value, cls, name):
    if not isinstance(value, cls):
        raise ContractError(f"{name} must be {cls.__name__}")


def _tuple(value, name):
    if type(value) is not tuple:
        raise ContractError(f"{name} must be an immutable tuple")


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode("utf-8")).hexdigest()


FEATURE_ORDER_SHA256 = _digest(NPB_FEATURES)


@dataclass(frozen=True)
class SourceFile:
    role: str
    path: str
    sha256: str | None

    def __post_init__(self):
        _text(self.role, "input role")
        _text(self.path, "input path")
        if not self.path.startswith("/") or any(p in ("", ".", "..") for p in self.path[1:].split("/")):
            raise ContractError("input path must be a literal absolute file path")
        if self.sha256 is not None:
            _sha(self.sha256, "input bytes")


@dataclass(frozen=True)
class SourceBinding:
    run: int
    segment: int
    sample: str
    inputs: tuple[SourceFile, ...]
    source_manifest_sha256: str
    retained_file_sha256: str
    reconstruction_sha256: str
    calibration_sha256: str
    input_identity_kind: str = "FILE_CONTENT_SHA256"
    input_tuple_record: str | None = None

    def __post_init__(self):
        _int(self.run, "run", 1)
        _int(self.segment, "segment")
        if self.sample not in SAMPLES:
            raise ContractError("unknown sample kind")
        _tuple(self.inputs, "input files")
        if not self.inputs:
            raise ContractError("original input files are required")
        for item in self.inputs:
            _typed(item, SourceFile, "input")
        if len({x.path for x in self.inputs}) != len(self.inputs) or len({x.role for x in self.inputs}) != len(self.inputs):
            raise ContractError("duplicate input path or role")
        if self.sample == "auau_data" and {x.role for x in self.inputs} != {"jets", "jetcalo"}:
            raise ContractError("AuAu DATA requires the exact jets/jetcalo pair")
        if self.input_identity_kind == 'FILE_CONTENT_SHA256':
            if any(item.sha256 is None for item in self.inputs) or self.input_tuple_record is not None:
                raise ContractError('file-content identity requires real checksums and no tuple substitute')
        elif self.input_identity_kind == 'FROZEN_TUPLE_RECORD_V1':
            # This binds paths/roles through the producer's sealed input record;
            # it explicitly does NOT claim a content checksum of each DST.
            record=self.input_tuple_record
            if (any(item.sha256 is not None for item in self.inputs)
                    or not isinstance(record,str) or len(record.encode('utf-8')) > 16384
                    or not record.endswith('\n') or record.count('\n') != 1
                    or '\r' in record or '\0' in record):
                raise ContractError('frozen tuple identity requires one literal record and null DST checksums')
            roles=(('jets','jetcalo') if self.sample.endswith('_data') else
                   ('calo_cluster','g4hits','jets','global','mbd_epd'))
            columns=record[:-1].split('\t')
            actual={item.role:item.path for item in self.inputs}
            if (len(columns)!=len(roles) or
                    actual!={role:path for role,path in zip(roles,columns) if path!='NONE'}
                    or not set(actual)<=set(roles)
                    or any(any(ch.isspace() for ch in path) for path in actual.values())):
                raise ContractError('frozen tuple paths or role order differ from source inputs')
            required={'jets','jetcalo'} if self.sample.endswith('_data') else {'jets','g4hits'}
            if not required<=set(actual):
                raise ContractError('frozen tuple is missing required source roles')
        else:
            raise ContractError('unknown input identity kind')
        for name in ("source_manifest_sha256", "retained_file_sha256", "reconstruction_sha256", "calibration_sha256"):
            _sha(getattr(self, name), name)

    @property
    def input_tuple_sha256(self):
        return (hashlib.sha256(self.input_tuple_record.encode('utf-8')).hexdigest()
                if self.input_identity_kind == 'FROZEN_TUPLE_RECORD_V1' else None)


@dataclass(frozen=True)
class EventKey:
    source: SourceBinding
    physical_event: int

    def __post_init__(self):
        _typed(self.source, SourceBinding, "source")
        _int(self.physical_event, "physical event", 0, 2**63 - 1)


@dataclass(frozen=True)
class CandidateKey:
    event: EventKey
    native_cluster_key: int

    def __post_init__(self):
        _typed(self.event, EventKey, "event key")
        _int(self.native_cluster_key, "rebuilt photon native key", 0, 2**32 - 1)


@dataclass(frozen=True)
class CandidateAnchor:
    cluster_et: float
    eta: float
    phi: float

    def __post_init__(self):
        for name in ("cluster_et", "eta", "phi"):
            _finite(getattr(self, name), name)
        if self.cluster_et < 0:
            raise ContractError("negative cluster ET")


@dataclass(frozen=True)
class SelectionFlags:
    reference_preselection_state: int
    active_preselection_state: int
    preselection_bitmask: int
    other_flags: tuple[tuple[str, int], ...] = ()

    def __post_init__(self):
        # These are opaque retained states, never inferred from the new score.
        for name in ("reference_preselection_state", "active_preselection_state"):
            if type(getattr(self, name)) is not int:
                raise ContractError("selection state must be a retained integer")
        _int(self.preselection_bitmask, "preselection bitmask", 0, 2**64 - 1)
        _tuple(self.other_flags, "other retained flags")
        names = set()
        for item in self.other_flags:
            _tuple(item, "flag item")
            if len(item) != 2:
                raise ContractError("flag item must be (name, integer)")
            name, value = item
            _text(name, "flag name")
            if name in names or name in ("reference_preselection_state", "active_preselection_state", "preselection_bitmask") or type(value) is not int:
                raise ContractError("duplicate, reserved or noninteger retained flag")
            names.add(name)


@dataclass(frozen=True)
class Measurement:
    state: str
    value: float | None = None
    reason: str = ""

    def __post_init__(self):
        if self.state not in ("VALID", "MISSING", "INVALID", "NOT_EVALUATED", "NOT_APPLICABLE"):
            raise ContractError("unknown measurement state")
        if self.value is not None:
            _finite(self.value, "measurement/raw sentinel")
        if self.state == "VALID":
            _finite(self.value, "valid measurement")
        else:
            _text(self.reason, "unavailable/invalid reason")
        if self.state == "MISSING" and self.value is not None:
            raise ContractError("missing data cannot contain an invented value")


@dataclass(frozen=True)
class ExistingField:
    """Finite retained value; unknown semantics are explicit and fail overlap.

    An adapter maps nonfinite native data to value=None without changing its base
    file. It must not assign this semantic hash based only on the branch name.
    """
    name: str
    value: float | None
    semantic_sha256: str | None

    def __post_init__(self):
        if self.name not in EVENT_FIELDS + CALO_FIELDS + CANDIDATE_FIELDS:
            raise ContractError("unknown augmentation field")
        if self.value is not None:
            _finite(self.value, "existing field")
        if self.semantic_sha256 is not None:
            _sha(self.semantic_sha256, "existing field semantics")


def _existing(fields, allowed):
    _tuple(fields, "existing fields")
    names = set()
    for item in fields:
        _typed(item, ExistingField, "existing field")
        if item.name not in allowed or item.name in names:
            raise ContractError("foreign or duplicate existing field")
        names.add(item.name)


@dataclass(frozen=True)
class RetainedEvent:
    key: EventKey
    existing: tuple[ExistingField, ...] = ()

    def __post_init__(self):
        _typed(self.key, EventKey, "retained event key")
        _existing(self.existing, EVENT_FIELDS + CALO_FIELDS)


@dataclass(frozen=True)
class RetainedCandidate:
    key: CandidateKey
    anchor: CandidateAnchor
    selection: SelectionFlags
    existing: tuple[ExistingField, ...] = ()

    def __post_init__(self):
        _typed(self.key, CandidateKey, "retained candidate key")
        _typed(self.anchor, CandidateAnchor, "retained ET/eta/phi")
        _typed(self.selection, SelectionFlags, "retained flags")
        _existing(self.existing, CANDIDATE_FIELDS)


@dataclass(frozen=True)
class TimingDefinition:
    definition: str
    config_sha256: str
    reconstruction_sha256: str
    calibration_sha256: str
    cluster_node: str
    tower_node: str
    sample_to_ns: float
    raw_units: str = "tower_time_samples"

    def __post_init__(self):
        if self.definition not in (NATIVE_TIMING, PPG12_TIMING):
            raise ContractError("unrecognized cluster timing population/definition")
        for name in ("config_sha256", "reconstruction_sha256", "calibration_sha256"):
            _sha(getattr(self, name), name)
        _text(self.cluster_node, "cluster node")
        _text(self.tower_node, "calibrated TowerInfo node")
        _finite(self.sample_to_ns, "sample-to-ns scale")
        if self.raw_units != "tower_time_samples" or self.sample_to_ns <= 0:
            raise ContractError("raw cluster timing must retain samples and an explicit positive scale")
        if self.definition == PPG12_TIMING and self.sample_to_ns != 17.6:
            raise ContractError("the audited PPG12 timing definition uses 17.6 ns/sample")

    @property
    def semantic_sha256(self):
        return _digest(asdict(self))

    @property
    def ownership(self):
        return "ALL_RAWCLUSTER_OWNED" if self.definition == NATIVE_TIMING else "RAWCLUSTER_OWNED_IN_7X7"

    @property
    def rejected_status_mask(self):
        return 0 if self.definition == NATIVE_TIMING else 1 << 5


@dataclass(frozen=True)
class TimingCapture:
    definition: TimingDefinition
    raw_time: Measurement
    denominator_energy: Measurement
    contributing_towers: int | None

    def __post_init__(self):
        _typed(self.definition, TimingDefinition, "timing definition")
        _typed(self.raw_time, Measurement, "raw timing")
        _typed(self.denominator_energy, Measurement, "timing denominator")
        if self.contributing_towers is not None:
            _int(self.contributing_towers, "contributing towers")
        if self.raw_time.state == "VALID":
            if self.denominator_energy.state != "VALID" or self.denominator_energy.value <= 0 or not self.contributing_towers:
                raise ContractError("valid timing requires a positive witnessed denominator and tower population")
            if self.definition.definition == NATIVE_TIMING and self.raw_time.value == -999:
                raise ContractError("native empty-denominator sentinel is not valid timing")
        elif self.raw_time.state == "NOT_APPLICABLE":
            raise ContractError("timing inapplicability is not established by this contract")
        if self.definition.definition == PPG12_TIMING and self.contributing_towers is not None and self.contributing_towers > 49:
            raise ContractError("PPG12 timing cannot use more than the 7x7 population")


@dataclass(frozen=True)
class MbdDefinition:
    reconstruction_sha256: str
    calibration_sha256: str
    config_sha256: str
    definition: str = "NATIVE_MBDOUT_T0_AND_ARM_TIMES_V1"
    units: str = "ns"
    arm_order: tuple[str, str] = ("south", "north")

    def __post_init__(self):
        for name in ("reconstruction_sha256", "calibration_sha256", "config_sha256"):
            _sha(getattr(self, name), name)
        if self.definition != "NATIVE_MBDOUT_T0_AND_ARM_TIMES_V1" or self.units != "ns" or self.arm_order != ("south", "north"):
            raise ContractError("MbdOut requires native ns values and arm0=south, arm1=north")

    @property
    def semantic_sha256(self):
        return _digest(asdict(self))


@dataclass(frozen=True)
class MbdCapture:
    definition: MbdDefinition
    native_is_valid: bool
    t0: Measurement
    south_time: Measurement
    north_time: Measurement
    south_npmt: int
    north_npmt: int

    def __post_init__(self):
        _typed(self.definition, MbdDefinition, "MBD definition")
        if type(self.native_is_valid) is not bool:
            raise ContractError("native MbdOut validity must be explicit")
        for value in (self.t0, self.south_time, self.north_time):
            _typed(value, Measurement, "native MBD time")
            if value.state == "NOT_APPLICABLE":
                raise ContractError("MbdOut inapplicability is not established by this contract")
            if value.state == "VALID" and value.value <= -9999:
                raise ContractError("MbdOut invalid sentinel is not a valid time")
        _int(self.south_npmt, "south npmt")
        _int(self.north_npmt, "north npmt")
        if self.native_is_valid != (self.t0.state == "VALID"):
            raise ContractError("native MbdOut validity contradicts t0 state")
        for hits, value in ((self.south_npmt, self.south_time), (self.north_npmt, self.north_time)):
            if hits == 0 and value.state == "VALID":
                raise ContractError("an empty MBD arm cannot be a valid measured arm time")


@dataclass(frozen=True)
class CaloDefinition:
    reconstruction_sha256: str
    calibration_sha256: str
    config_sha256: str
    require_isgood: bool = True
    definition: str = "RJ_CALIB_SIGNED_FINITE_DOUBLE_SUM_FLOAT_STORE_V1"
    units: str = "GeV"
    tower_nodes: tuple[str, ...] = ("TOWERINFO_CALIB_CEMC", "TOWERINFO_CALIB_HCALIN", "TOWERINFO_CALIB_HCALOUT")

    def __post_init__(self):
        for name in ("reconstruction_sha256", "calibration_sha256", "config_sha256"):
            _sha(getattr(self, name), name)
        if (type(self.require_isgood) is not bool or self.units != "GeV" or
                self.definition != "RJ_CALIB_SIGNED_FINITE_DOUBLE_SUM_FLOAT_STORE_V1" or
                self.tower_nodes != ("TOWERINFO_CALIB_CEMC", "TOWERINFO_CALIB_HCALIN", "TOWERINFO_CALIB_HCALOUT")):
            raise ContractError("calorimeter definition must bind native calibrated nodes, precision and isGood policy")

    @property
    def semantic_sha256(self):
        return _digest(asdict(self))


@dataclass(frozen=True)
class CaloCapture:
    definition: CaloDefinition
    available_mask: int
    valid_mask: int
    values: tuple[Measurement, ...]

    def __post_init__(self):
        _typed(self.definition, CaloDefinition, "calorimeter definition")
        _int(self.available_mask, "calorimeter available mask", 0, 7)
        _int(self.valid_mask, "calorimeter valid mask", 0, 7)
        if self.valid_mask & ~self.available_mask:
            raise ContractError("calorimeter validity contradicts container availability")
        _tuple(self.values, "calorimeter values")
        if len(self.values) != 4:
            raise ContractError("calorimeter capture needs three layers and a total")
        for index, value in enumerate(self.values):
            _typed(value, Measurement, "calorimeter value")
            if value.state not in ("VALID", "MISSING", "INVALID"):
                raise ContractError("calorimeter value cannot be inapplicable or unevaluated")
            available = bool(self.available_mask & (1 << index)) if index < 3 else self.available_mask == 7
            valid = bool(self.valid_mask & (1 << index)) if index < 3 else self.valid_mask == 7
            if not available and (value.state != "MISSING" or value.value is not None):
                raise ContractError("missing calorimeter container cannot supply a scalar")
            if index < 3 and valid != (value.state == "VALID"):
                raise ContractError("calorimeter valid mask contradicts layer measurement")
            if index == 3 and value.state == "VALID" and not valid:
                raise ContractError("calorimeter total needs three valid layers")
            if index == 3 and valid and value.value is not None and value.state != "VALID":
                raise ContractError("finite calorimeter total with three valid layers must be valid")


def _float32(value):
    try:
        return struct.unpack("!f", struct.pack("!f", value))[0]
    except (OverflowError, struct.error) as exc:
        raise ContractError("NPB float32 domain overflow") from exc


def _float32_bits_value(bits):
    if not isinstance(bits, str) or not re.fullmatch(r"[0-9a-f]{8}", bits):
        raise ContractError("native float32 bits must be eight lowercase hex digits")
    return struct.unpack("!f", bytes.fromhex(bits))[0]


@dataclass(frozen=True)
class NativeFloat32Bound:
    """Exact native bound; NaN and either infinity independently disable it.

    Hex encodes the IEEE-754 binary32 word in big-endian order. No nonfinite
    JSON number is introduced, and DISABLED is distinct from unknown/missing.
    The definition's configuration/evidence hashes must bind these actual bits.
    """
    state: str
    bits: str

    def __post_init__(self):
        value = _float32_bits_value(self.bits)
        if self.state not in ("FINITE", "DISABLED_NONFINITE"):
            raise ContractError("unknown native NPB bound state")
        if math.isfinite(value) != (self.state == "FINITE"):
            raise ContractError("native NPB bound state contradicts captured float bits")

    @property
    def value(self):
        return _float32_bits_value(self.bits)


@dataclass(frozen=True)
class ScoringObjectWitness:
    """The object actually passed to the NPB accessor, not an inferred alias.

    SAME_OBJECT_V1 requires exact retained identity/anchor. The native fallback
    accepts the first kinematic match, not necessarily exact ET/eta equality.
    Its receipt must prove the named scoring node, selected key and encounter
    ordinal, including absence of an earlier passing row. Hash syntax alone is
    not that proof. Input bits come from the actual float accessor, not a Python
    recomputation of energy/cosh or an arbitrary anchor-equality tolerance.
    """
    key: CandidateKey
    cluster_node: str
    anchor: CandidateAnchor
    input_et_bits: str
    input_eta_bits: str
    binding_policy: str
    evidence_sha256: str
    matched_encounter_ordinal: int | None = None

    def __post_init__(self):
        _typed(self.key, CandidateKey, "NPB scoring-object key")
        _typed(self.anchor, CandidateAnchor, "NPB scoring-object stored anchor")
        _text(self.cluster_node, "NPB scoring node")
        _sha(self.evidence_sha256, "NPB scoring-object binding evidence")
        et = _float32_bits_value(self.input_et_bits)
        eta = _float32_bits_value(self.input_eta_bits)
        if not math.isfinite(et) or et < 0 or not math.isfinite(eta):
            raise ContractError("NPB scoring kinematic input bits must be finite")
        if eta != self.anchor.eta:
            raise ContractError("NPB eta accessor bits differ from scoring object's stored eta")
        if self.binding_policy == "SAME_OBJECT_V1":
            if self.matched_encounter_ordinal is not None:
                raise ContractError("same-object NPB is not a fallback encounter match")
        elif self.binding_policy == "NATIVE_FIRST_KINEMATIC_MATCH_V1":
            _int(self.matched_encounter_ordinal, "first matched scoring-object ordinal", 0, 2**63 - 1)
        else:
            raise ContractError("unknown NPB scoring-object binding policy")

    @property
    def input_anchor(self):
        return CandidateAnchor(_float32_bits_value(self.input_et_bits),
                               _float32_bits_value(self.input_eta_bits), self.anchor.phi)

    def validate_binding(self, retained_key, retained_anchor, definition):
        if self.key.event != retained_key.event or self.cluster_node != definition.cluster_node:
            raise ContractError("foreign NPB scoring source/event/node")
        if self.binding_policy == "SAME_OBJECT_V1":
            if self.key != retained_key or self.anchor != retained_anchor:
                raise ContractError("same-object NPB identity/anchor differs from retained candidate")
            return
        # Literal native findMatchedPhotonByKinematics predicate (not a new
        # tolerance). First-match/order validity is an external receipt gate.
        dphi = self.anchor.phi - retained_anchor.phi
        # Native stored angles must be canonical; avoid normalizing invented
        # arbitrary-angle inputs with a long-running loop.
        # PhotonClusterv1 stores floats: binary32 pi rounds just outside the
        # mathematical double-precision interval and is still a native angle.
        phi_limit = _float32(math.pi)
        if not (-phi_limit <= self.anchor.phi <= phi_limit and
                -phi_limit <= retained_anchor.phi <= phi_limit):
            raise ContractError("NPB fallback needs canonical native phi witnesses")
        if dphi >= math.pi:
            dphi -= 2 * math.pi
        if dphi < -math.pi:
            dphi += 2 * math.pi
        if not (abs(self.anchor.eta - retained_anchor.eta) < 1e-5 and
                abs(dphi) < 1e-5 and
                abs(self.anchor.cluster_et - retained_anchor.cluster_et) < 1e-4):
            raise ContractError("NPB fallback contradicts native first-match kinematic predicate")


@dataclass(frozen=True)
class NPBDefinition:
    config_sha256: str
    applicability_sha256: str
    feature_semantics_sha256: str
    reconstruction_sha256: str
    calibration_sha256: str
    cluster_node: str
    applicability_name: str
    samples: tuple[str, ...]
    et_min: float | None
    et_max: float | None
    abs_eta_max: float | None
    enabled: bool
    score_field: str = "npb_score"
    model_name: str = "npb_score_split"
    model_sha256: str = SPLIT_MODEL_SHA256
    training_config_sha256: str = SPLIT_TRAINING_CONFIG_SHA256
    feature_names: tuple[str, ...] = NPB_FEATURES
    feature_order_sha256: str = FEATURE_ORDER_SHA256
    domain_policy: str = "EXPLICIT_CLOSED_V1"
    configured_cut: float | None = None
    native_bounds: tuple[NativeFloat32Bound, ...] = ()

    def __post_init__(self):
        for name in ("config_sha256", "applicability_sha256", "feature_semantics_sha256", "reconstruction_sha256", "calibration_sha256"):
            _sha(getattr(self, name), name)
        for name in ("cluster_node", "applicability_name"):
            _text(getattr(self, name), name)
        if (self.model_name != "npb_score_split" or self.model_sha256 != SPLIT_MODEL_SHA256
                or self.training_config_sha256 != SPLIT_TRAINING_CONFIG_SHA256):
            raise ContractError("this adapter is bound to the audited split NPB model/training config")
        _tuple(self.feature_names, "NPB feature names")
        if self.feature_names != NPB_FEATURES or self.feature_order_sha256 != FEATURE_ORDER_SHA256:
            raise ContractError("NPB requires the exact 25-input order/hash")
        _tuple(self.samples, "NPB sample applicability")
        if not self.samples or len(set(self.samples)) != len(self.samples) or any(x not in SAMPLES for x in self.samples):
            raise ContractError("invalid explicit NPB sample applicability")
        if self.score_field not in ("npb_score", "auau_npb_score"):
            raise ContractError("NPB score field must be named explicitly")
        if self.domain_policy not in ("EXPLICIT_CLOSED_V1", "NAMED_FLOAT32_V1", "PRIMARY_UNGATED_V1"):
            raise ContractError("unknown NPB producer domain policy")
        if type(self.enabled) is not bool:
            raise ContractError("NPB enablement must be explicit")
        _tuple(self.native_bounds, "native NPB bounds")
        bound_names = ("et_min", "et_max", "abs_eta_max")
        if self.native_bounds:
            if self.domain_policy != "NAMED_FLOAT32_V1" or len(self.native_bounds) != 3:
                raise ContractError("native NPB bounds require the three named-float32 gates")
            for name, bound in zip(bound_names, self.native_bounds):
                _typed(bound, NativeFloat32Bound, "native NPB bound")
                value = getattr(self, name)
                if bound.state == "DISABLED_NONFINITE":
                    if value is not None:
                        raise ContractError("disabled native NPB bound must not invent a finite cutoff")
                else:
                    _finite(value, name)
                    if _float32(value) != bound.value:
                        raise ContractError("native NPB bound bits contradict the configured value")
        else:
            for name in bound_names:
                value = getattr(self, name)
                if value is None and self.domain_policy == "PRIMARY_UNGATED_V1":
                    continue  # Explicitly unused by the primary evaluator.
                _finite(value, name)
        if self.domain_policy == "EXPLICIT_CLOSED_V1" and (
                self.et_min < 0 or self.et_max < self.et_min or self.abs_eta_max < 0):
            raise ContractError("invalid explicit NPB kinematic applicability")
        if self.configured_cut is not None:
            _finite(self.configured_cut, "NPB configured strict-greater cut")
        # JSON emitters may serialize 6.0 as 6. These are identical physical
        # bounds and native float32 gates; canonicalize their Python types so
        # a roundtrip cannot change the semantic hash or break finite overlap.
        for name in (*bound_names, "configured_cut"):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, float(getattr(self, name)))

    @property
    def semantic_sha256(self):
        return _digest(asdict(self))

    def contains(self, sample, anchor):
        if not self.enabled or sample not in self.samples:
            return False
        if self.domain_policy == "PRIMARY_UNGATED_V1":
            return True  # Audited primary evaluator has no additional domain gate.
        if self.domain_policy == "NAMED_FLOAT32_V1":
            et, eta = map(_float32, (anchor.cluster_et, anchor.eta))
            if self.native_bounds:
                low, high, eta_max = (x.value for x in self.native_bounds)
            else:
                low, high, eta_max = map(_float32, (self.et_min, self.et_max, self.abs_eta_max))
            return ((not math.isfinite(low) or et > low or
                     abs(_float32(et - low)) < _float32(1e-6)) and
                    (not math.isfinite(high) or et <= high) and
                    (not math.isfinite(eta_max) or abs(eta) <= eta_max))
        return self.et_min <= anchor.cluster_et <= self.et_max and abs(anchor.eta) <= self.abs_eta_max


def _allowed_origins(name):
    if name == "vertexz":
        return {"NATIVE_PRODUCER", "RETAINED_VERTEX_EQUIVALENCE_VERIFIED"}
    if name == "e11_over_e22" or name.startswith("e22_over_"):
        return {"NATIVE_PRODUCER", "REPLAYED_RAWCLUSTER_MAP_VERIFIED"}
    if name in ("e11_over_e33", "e32_over_e35"):
        return {"NATIVE_PRODUCER", "NPB_FLOAT_RATIO_VERIFIED"}
    if "_over_" in name or name in ("cluster_w32", "cluster_w52", "cluster_w72"):
        return {"NATIVE_PRODUCER", "REPLAYED_CALIBRATED_GRID_VERIFIED"}
    return {"NATIVE_PRODUCER", "RETAINED_NATIVE_VERIFIED"}


@dataclass(frozen=True)
class FeatureInput:
    name: str
    measurement: Measurement
    origin: str
    evidence_sha256: str | None

    def __post_init__(self):
        if self.name not in NPB_FEATURES:
            raise ContractError("unknown NPB input")
        _typed(self.measurement, Measurement, "NPB feature value")
        if self.origin == "UNAVAILABLE":
            if self.measurement.state == "VALID" or self.measurement.value is not None or self.evidence_sha256 is not None:
                raise ContractError("unavailable NPB input cannot be zero-filled or declared valid")
        else:
            if self.origin not in _allowed_origins(self.name):
                raise ContractError(f"unproven NPB origin/semantics for {self.name}")
            _sha(self.evidence_sha256, "NPB input capture/replay/equivalence evidence")


@dataclass(frozen=True)
class NPBCapture:
    definition: NPBDefinition
    inputs: tuple[FeatureInput, ...]
    applicability_state: str
    evaluation_state: str
    score: Measurement
    scoring_object: ScoringObjectWitness | None = None

    def __post_init__(self):
        _typed(self.definition, NPBDefinition, "NPB definition")
        _tuple(self.inputs, "NPB inputs")
        for item in self.inputs:
            _typed(item, FeatureInput, "NPB input")
        if tuple(x.name for x in self.inputs) != NPB_FEATURES:
            raise ContractError("NPB input witnesses must contain exactly the ordered 25 names")
        _typed(self.score, Measurement, "NPB score")
        if self.scoring_object is not None:
            _typed(self.scoring_object, ScoringObjectWitness, "NPB scoring-object witness")
        if self.applicability_state not in ("APPLICABLE", "NOT_APPLICABLE", "UNVERIFIED"):
            raise ContractError("unknown NPB applicability state")
        if self.evaluation_state not in ("EVALUATED", "EVALUATED_INVALID", "NOT_EVALUATED", "INPUTS_UNAVAILABLE", "NOT_APPLICABLE"):
            raise ContractError("unknown NPB evaluation state")
        if self.evaluation_state == "EVALUATED":
            if (self.applicability_state != "APPLICABLE" or self.score.state != "VALID"
                    or any(x.measurement.state != "VALID" for x in self.inputs)
                    or not 0 <= self.score.value <= 1):
                raise ContractError("evaluated NPB requires applicability, all25 valid inputs and a nonsentinel probability")
        elif self.evaluation_state == "EVALUATED_INVALID":
            if self.applicability_state != "APPLICABLE" or self.score.state != "INVALID":
                raise ContractError("invalid evaluated NPB must retain applicable/invalid state")
            if (all(x.measurement.state == "VALID" for x in self.inputs) and
                    self.score.value is not None and 0 <= self.score.value <= 1):
                raise ContractError("invalid evaluated NPB needs invalid inputs or invalid output")
        elif self.evaluation_state == "NOT_APPLICABLE":
            if self.applicability_state != "NOT_APPLICABLE" or self.score.state != "NOT_APPLICABLE":
                raise ContractError("NPB inapplicability states disagree")
        else:
            if self.applicability_state == "NOT_APPLICABLE" or self.score.state not in ("NOT_EVALUATED", "MISSING", "INVALID"):
                raise ContractError("unevaluated NPB cannot claim a valid score")
            if self.evaluation_state == "INPUTS_UNAVAILABLE" and all(x.measurement.state == "VALID" for x in self.inputs):
                raise ContractError("INPUTS_UNAVAILABLE contradicts all25 valid inputs")
        if self.evaluation_state not in ("EVALUATED", "EVALUATED_INVALID") and self.score.value not in (None, -1):
            raise ContractError("unevaluated NPB must retain only the native -1 sentinel or absence")

    def validate_candidate(self, key, anchor):
        if self.scoring_object is not None:
            self.scoring_object.validate_binding(key, anchor, self.definition)
            input_anchor = self.scoring_object.input_anchor
        else:
            if self.definition.domain_policy != "EXPLICIT_CLOSED_V1":
                raise ContractError("native NPB policy needs an explicit scoring-object witness")
            input_anchor = anchor  # Legacy explicit exact-equality contract only.
        inside = self.definition.contains(key.event.source.sample, input_anchor)
        if self.applicability_state != "UNVERIFIED" and inside != (self.applicability_state == "APPLICABLE"):
            raise ContractError("NPB applicability contradicts explicit producer domain")
        # No event-vertex equality is inferred here; vertex origin needs evidence.
        for item, expected in zip(self.inputs[:2], (input_anchor.cluster_et, input_anchor.eta)):
            if item.measurement.state == "VALID" and item.measurement.value != expected:
                raise ContractError("NPB ET/eta inputs contradict bound scoring-object accessor")

    @property
    def pass_measurement(self):
        """New explicit discriminator only; never replace retained selection flags."""
        if self.evaluation_state == "NOT_APPLICABLE":
            return Measurement("NOT_APPLICABLE", reason="outside the explicitly bound NPB domain")
        if self.evaluation_state == "EVALUATED_INVALID":
            return Measurement("INVALID", reason="native Compute ran but its inputs or output were invalid")
        if self.evaluation_state != "EVALUATED":
            return Measurement("NOT_EVALUATED", reason="no valid evaluated NPB score")
        if self.definition.configured_cut is None:
            return Measurement("MISSING", reason="configured NPB cut was not captured/bound")
        return Measurement("VALID", int(self.score.value > self.definition.configured_cut))


@dataclass(frozen=True)
class Missing:
    reason: str

    def __post_init__(self):
        _text(self.reason, "missing donor section reason")


@dataclass(frozen=True)
class EventDonor:
    key: EventKey
    capture_sha256: str
    mbd: MbdCapture | Missing
    calo: CaloCapture | Missing = Missing("calorimeter donor not requested/captured")

    def __post_init__(self):
        _typed(self.key, EventKey, "donor event key")
        _sha(self.capture_sha256, "donor capture receipt")
        if not isinstance(self.mbd, (MbdCapture, Missing)):
            raise ContractError("MBD donor must be captured or explicitly missing")
        if not isinstance(self.calo, (CaloCapture, Missing)):
            raise ContractError("calorimeter donor must be captured or explicitly missing")


@dataclass(frozen=True)
class CandidateDonor:
    key: CandidateKey
    anchor: CandidateAnchor
    selection: SelectionFlags
    capture_sha256: str
    timing: TimingCapture | Missing
    npb: NPBCapture | Missing

    def __post_init__(self):
        _typed(self.key, CandidateKey, "donor candidate key")
        _typed(self.anchor, CandidateAnchor, "donor ET/eta/phi")
        _typed(self.selection, SelectionFlags, "donor flags")
        _sha(self.capture_sha256, "donor capture receipt")
        if not isinstance(self.timing, (TimingCapture, Missing)) or not isinstance(self.npb, (NPBCapture, Missing)):
            raise ContractError("donor sections must be captured or explicitly missing")


@dataclass(frozen=True)
class OverlayField:
    name: str
    measurement: Measurement
    semantic_sha256: str


@dataclass(frozen=True)
class OverlayRow:
    retained: RetainedEvent | RetainedCandidate
    donor: EventDonor | CandidateDonor | None
    fields: tuple[OverlayField, ...]

    @property
    def complete(self):
        return self.donor is not None and all(x.measurement.state in ("VALID", "NOT_APPLICABLE") for x in self.fields)


@dataclass(frozen=True)
class OverlayReport:
    source: SourceBinding
    events: tuple[OverlayRow, ...]
    candidates: tuple[OverlayRow, ...]
    missing_event_keys: tuple[EventKey, ...]
    missing_candidate_keys: tuple[CandidateKey, ...]
    scientific_acceptance: bool = False
    capture_source: SourceBinding | None = None

    @property
    def complete(self):
        return bool(self.events) and all(x.complete for x in self.events + self.candidates)

    @property
    def unresolved(self):
        return tuple((row.retained.key, f.name, f.measurement.state, f.measurement.reason)
                     for row in self.events + self.candidates for f in row.fields
                     if f.measurement.state not in ("VALID", "NOT_APPLICABLE"))


def _index(rows, cls, label):
    result = {}
    for row in rows:
        _typed(row, cls, label)
        if row.key in result:
            raise ContractError(f"duplicate {label} key")
        result[row.key] = row
    return result


def _overlap(retained, fields):
    original = {x.name: x for x in retained.existing}
    for field in fields:
        old = original.get(field.name)
        if old is None or old.value is None or field.measurement.value is None:
            continue
        # Even finite native sentinels stay immutable. A new semantic definition
        # belongs in a different explicitly named contract, never a branch alias.
        if old.semantic_sha256 != field.semantic_sha256:
            raise ContractError(f"unbound/different finite-overlap semantics: {field.name}")
        if old.value != field.measurement.value:
            raise ContractError(f"contradictory finite retained/donor value: {field.name}")


def _field(name, value, definition):
    return OverlayField(name, value, definition.semantic_sha256)


def reconcile(source, events, candidates, event_donors, candidate_donors, *, timing, mbd, npb, calo=None,
              capture_source=None):
    """Produce immutable overlays, or fail on ambiguity/contradiction.

    Timing, MBD and NPB are requested; calorimeter augmentation additionally
    requires an explicit CaloDefinition. Missing requested sections/rows remain explicit and
    report.complete is false. NOT_APPLICABLE is accepted only for NPB after
    checking the hash-bound producer domain; it never changes native selection.
    Empty candidate populations are permitted, but an empty event scope is not.
    No entry ordinal participates in any join. Exact ET/eta/phi equality is
    required; callers cannot weaken it with a tolerance in this adapter.
    """
    _typed(source, SourceBinding, "scope source")
    definition_source = source if capture_source is None else capture_source
    _typed(definition_source, SourceBinding, "capture source")
    # A changed producer may add captures; it cannot silently change input
    # data or the calibration basis. Original retained identities stay scoped
    # to source, while added field definitions name the actual capture build.
    if any(getattr(source, name) != getattr(definition_source, name) for name in
           ('run', 'segment', 'sample', 'inputs', 'calibration_sha256')):
        raise ContractError('capture source changes input identity or calibration basis')
    definitions = [(timing, TimingDefinition), (mbd, MbdDefinition), (npb, NPBDefinition)]
    if calo is not None:
        definitions.append((calo, CaloDefinition))
    for definition, cls in definitions:
        _typed(definition, cls, "expected definition")
        if (definition.reconstruction_sha256 != definition_source.reconstruction_sha256
                or definition.calibration_sha256 != definition_source.calibration_sha256):
            raise ContractError("expected definition is not bound to source reconstruction/calibration")
    base_events = _index(events, RetainedEvent, "retained event")
    base_candidates = _index(candidates, RetainedCandidate, "retained candidate")
    donors_e = _index(event_donors, EventDonor, "event donor")
    donors_c = _index(candidate_donors, CandidateDonor, "candidate donor")
    if not base_events:
        raise ContractError("empty retained event scope")
    if any(key.source != source for key in base_events):
        raise ContractError("foreign retained source/hash")
    if any(key.event not in base_events for key in base_candidates):
        raise ContractError("candidate has no retained physical event/source")
    if not donors_e.keys() <= base_events.keys() or not donors_c.keys() <= base_candidates.keys():
        raise ContractError("foreign donor source/hash/event/object")
    out_e, out_c = [], []
    for key, base in sorted(base_events.items(), key=lambda x: x[0].physical_event):
        donor = donors_e.get(key)
        captured = donor.mbd if donor else Missing("no donor for exact source/physical event key")
        if isinstance(captured, Missing):
            values = (Measurement("MISSING", reason=captured.reason),) * 3
        else:
            if captured.definition != mbd:
                raise ContractError("MbdOut donor definition/config differs from requested semantics")
            values = captured.t0, captured.south_time, captured.north_time
        fields = tuple(_field(name, value, mbd) for name, value in zip(EVENT_FIELDS, values))
        if calo is not None:
            capture = donor.calo if donor else Missing("no calorimeter donor for exact source/physical event key")
            if isinstance(capture, Missing):
                calo_values = (Measurement("MISSING", reason=capture.reason),) * 4
            else:
                if capture.definition != calo:
                    raise ContractError("calorimeter donor definition/config differs from requested semantics")
                calo_values = capture.values
            fields += tuple(_field(name, value, calo) for name, value in zip(CALO_FIELDS, calo_values))
        _overlap(base, fields)
        out_e.append(OverlayRow(base, donor, fields))
    for key, base in sorted(base_candidates.items(), key=lambda x: (x[0].event.physical_event, x[0].native_cluster_key)):
        donor = donors_c.get(key)
        if donor and (donor.anchor != base.anchor or donor.selection != base.selection):
            raise ContractError("donor rebuilt candidate/ET/eta/phi or immutable selection flags differ")
        time_capture = donor.timing if donor else Missing("no donor for exact source/event/rebuilt photon key")
        npb_capture = donor.npb if donor else Missing("no donor for exact source/event/rebuilt photon key")
        if isinstance(time_capture, Missing):
            time_value = Measurement("MISSING", reason=time_capture.reason)
        else:
            if time_capture.definition != timing:
                raise ContractError("cluster timing donor population/config/units differ from requested semantics")
            time_value = time_capture.raw_time
        if isinstance(npb_capture, Missing):
            score = Measurement("MISSING", reason=npb_capture.reason)
        else:
            if npb_capture.definition != npb:
                raise ContractError("NPB donor model/config/applicability/feature semantics differ")
            npb_capture.validate_candidate(key, base.anchor)
            score = npb_capture.score
        fields = (_field("cluster_time_raw", time_value, timing), _field(npb.score_field, score, npb))
        _overlap(base, fields)
        out_c.append(OverlayRow(base, donor, fields))
    return OverlayReport(source, tuple(out_e), tuple(out_c),
                         tuple(row.retained.key for row in out_e if row.donor is None),
                         tuple(row.retained.key for row in out_c if row.donor is None),
                         capture_source=capture_source)
