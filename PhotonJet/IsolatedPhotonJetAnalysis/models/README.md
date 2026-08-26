# Models

The two compact H70 TMVA files in `models/artifacts/` are redistribution-safe,
hash-bound references. Their `PhotonJetModelReferenceV1` records bind only the
artifact bytes, collision system, feature-contract pointer, feature count, and
score direction.

They are not consumed by the current tree-starting V1 examples or CLI. Public
model inference, score reproduction, training-dataset/split provenance,
held-out validation, working-point derivation, and production equivalence are
all `NOT_RUN`; the reference records must not be interpreted as model cards or
physics validation.
