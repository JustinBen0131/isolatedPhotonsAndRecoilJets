# Producer integration

The current production implementation remains the golden oracle while its
configuration and system-specific code are separated behind differential
tests. It is not copied here under a new class name before that proof exists.

This draft starts at the accepted public tree boundary. The private migration
adapter used to establish that boundary is intentionally absent. The eventual
source files in this directory are:

- `PhotonJetAnalysis.{h,cc}`: compatibility facade and event loop;
- `PhotonJetConfig.{h,cc}`: strict typed configuration;
- `PhotonJetIdentity.{h,cc}`: occurrence-safe identities;
- `PhotonJetRecords.{h,cc}`: public event records;
- `PhotonJetTreeWriter.{h,cc}`: fail-closed tree writer;
- `PhotonJetSelections.{h,cc}`: shared photon and recoil definitions;
- `PhotonJetAuAuAdapter.{h,cc}`: heavy-ion nodes, trigger, event gate,
  centrality, and background state.

No mini-DST producer or direct-reference equivalence claim is made by the current draft. See
`docs/release_status.md`.
