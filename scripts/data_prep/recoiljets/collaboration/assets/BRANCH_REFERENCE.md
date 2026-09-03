# PhotonJetTrees branch reference

All numeric branches use standard ROOT integer or floating-point types.
Object identities are stored as `_hi` and `_lo` unsigned 64-bit components.
Within one part, `source_file_index`, `event_id_hi`, and `event_id_lo` identify
an event. When using a `TChain`, include `TChain::GetTreeNumber()` as the file
component of an event or object key. The supplied examples do this.

## `events`

One entry per event.

- identity: `source_file_index`, `source_entry`, `event_id_hi`, `event_id_lo`;
- event record: `run`, `event_sequence`, `physical_event_sequence`;
- trigger: `trigger_bits`, `live_trigger_bits`, `scaled_trigger_bits`,
  `scaled_bit30`; Au+Au DATA additionally stores `scaled_bit22`,
  `minimum_bias_pass`, and `nominal_event_selection_pass`;
- analysis values: `vertex_z`, `centrality`, `event_weight`,
  `total_calo_energy`, `terminal_status`.

`event_weight` is copied verbatim from the producer input. It is not guaranteed
to contain every normalization factor required by a downstream analysis. For
Au+Au embedded-inclusive Jet12/20/30/40, it contains neither the accepted
source-slice stitch factor (`ownership_effective_cross_section / Npass` in the
matching half-open maximum-truth-jet window) nor the separately applied
centrality factor. The package has no global per-row source-slice sidecar, so a
complete nominal Au+Au embedded-inclusive weight cannot be reconstructed from
these v1 Trees alone. Do not treat `source_file_index`, part order, filenames,
or `event_weight` as a substitute. A nominal consumer must use a separately
supplied hash-bound composite weight payload or report
`NOT_READY_FROM_PACKAGE_ALONE`.

`scaled_bit30` is meaningful for the p+p data trigger selection. Simulation
and samples without that trigger record use the documented unavailable value.

For Au+Au DATA, `run` is the positive run number recovered from the exact
source-bound event-gate witness. `scaled_trigger_bits` is the recorded GL1
scaled vector, `scaled_bit22` is its direct bit-22 value, and
`minimum_bias_pass` is the recorded `MinimumBiasInfo::isAuAuMinimumBias()`
decision. Every retained Au+Au DATA row has `scaled_bit22 == 1`,
`minimum_bias_pass == 1`, and `nominal_event_selection_pass == 1`; the same
event-level values are repeated consistently in all five run-bearing trees.
See `validation/AUAU_RUN_GATE_PROVENANCE.json` for the release-bound source,
gate, and aggregate receipts.

## `eventTree`

One entry per event. It contains every `events` branch plus arrays that can be
looped without joining separate TTrees.

Photon arrays:

`photon_candidate_id_hi`, `photon_candidate_id_lo`, `photon_encounter_ordinal`, `photon_et`, `photon_eta`,
`photon_phi`, `photon_bdt_score`, `photon_bdt_tight_threshold`,
`photon_bdt_nontight_low_threshold`, `photon_bdt_nontight_high_threshold`,
`photon_bdt_is_tight`, `photon_bdt_is_nontight`,
`photon_bdt_is_not_tight`, `photon_iso_r03`, `photon_iso_r03_threshold`,
`photon_iso_r03_nonisolated_threshold`, `photon_iso_r04`,
`photon_iso_r04_threshold`, `photon_iso_r04_nonisolated_threshold`.

Jet arrays:

`jet_id_hi`, `jet_id_lo`, `jet_pt`, `jet_raw_pt`, `jet_eta`, `jet_phi`,
`jet_radius`, `jet_mass`, `jet_area`.

Pair arrays:

`pair_photon_index`, `pair_jet_index`, `pair_delta_phi`, `pair_xjgamma`,
`pair_recoil_state`. The two index arrays point into the photon and jet arrays
in the same event entry.

Counts are stored as `nphotons`, `njets`, and `npairs`.

The convenience branches `leader_A_r04_index` through `leader_D_r04_index`
contain event-local indices for the established tight/bounded-non-tight score
definitions and R=0.4 isolation. `leader_C_r04_complement_index` and
`leader_D_r04_complement_index` instead use the full tight-score complement.
These convenience indices are computed over the complete stored photon list.
Recompute leaders from the photon arrays after applying the analysis photon
pT, eta, vertex, and centrality selection, or whenever using another isolation
or score definition. This preserves independent event-leading choices in A,
B, C, and D without making a hidden acceptance choice.

## `photons`

One entry per photon candidate with complete, finite inputs for the included
model. It repeats the `events` branches and adds:

- identity and kinematics: `candidate_id_hi`, `candidate_id_lo`,
  `photon_encounter_ordinal`, `photon_et`, `photon_eta`, `photon_phi`;
- score and established selections: `bdt_score`, `bdt_tight_threshold`,
  `bdt_nontight_low_threshold`, `bdt_nontight_high_threshold`,
  `bdt_is_tight`, `bdt_is_nontight`, `bdt_is_not_tight`,
  `bdt_input_count`, `bdt_input_00` through `bdt_input_13`;
- R=0.3 isolation: `iso_r03`, `iso_r03_threshold`,
  `iso_r03_nonisolated_threshold`, `iso_r03_pass`;
- R=0.4 isolation: `iso_r04`, `iso_r04_threshold`,
  `iso_r04_nonisolated_threshold`, `iso_r04_pass`;
- simulation match summary: `truth_matched`, `truth_barcode`;
- selected native shower measurements: `native_weta_cogx`,
  `native_wphi_cogx`, `native_weta33_cogx`, `native_wphi33_cogx`,
  `native_e11_over_e33`, `native_e32_over_e35`, `native_et1`.

The p+p score uses 11 ordered inputs and the Au+Au score uses 14. Unused
input slots are not part of that model's evaluation.

## `jets`

One entry per reconstructed anti-kT R=0.4 jet. It repeats the `events`
branches and adds `jet_id_hi`, `jet_id_lo`, `jet_radius`, `jet_raw_pt`,
`jet_pt`, `jet_eta`, `jet_phi`, `jet_mass`, `jet_area`,
`jet_quality_bitmask`, and `jet_order`.

## `photonJets`

One entry per stored photon-jet pair whose photon is present in `photons`. It
repeats all `events`, `photons`, and `jets` branches and adds:

`pair_id_hi`, `pair_id_lo`, `photon_index`, `jet_index`, `delta_phi`,
`xjgamma`, `recoil_state`, `photon_rank`, `jet_rank`,
`wrong_photon_class`, `wrong_recoil_class`.

`photon_index` and `jet_index` are local indices within the corresponding
event. `xjgamma` is the stored `jet_pt / photon_et` pair observable.
`delta_phi` is the absolute wrapped photon-jet azimuthal separation,
`abs(atan2(sin(photon_phi - jet_phi), cos(photon_phi - jet_phi)))`, in
`[0, pi]`. Both relations are enforced in the source and release validators.

## Simulation truth trees

`truthPhotons` contains `source_file_index`, event and truth-photon identities,
`truth_photon_pt`, `truth_photon_eta`, `truth_photon_phi`, `prompt_class`,
`source_role`, `generator_barcode`, and `truth_isolation`.

`truthJets` contains `source_file_index`, event and truth-jet identities,
`truth_jet_radius`, `truth_jet_pt`, `truth_jet_eta`, and `truth_jet_phi`.

`recoTruthLinks` contains `source_file_index`, event and link identities,
reconstructed and truth object identities and types, event-local `reco_index`
and `truth_index`, `match_metric`, and `link_class`. A local index is `-1` when
that side of the link is absent. Analyses must state any link-class or
best-match requirement rather than inferring one from file order.
