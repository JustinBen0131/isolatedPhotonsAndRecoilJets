# Quick start

Use the same ROOT environment that provides `uproot`, `awkward`, and the
sPHENIX model runtime:

```bash
python3 -m pip install -e '.[plots]' -c requirements-constraints.txt

photonjet trees validate \
  --input PhotonJetTrees_v1.root \
  --model-input-count 11

photonjet histogram \
  --input PhotonJetTrees_v1.root \
  --region A \
  --dataset config/datasets/pp_engineering_fixture.yaml \
  --output xjgamma.json

photonjet plot render \
  --histogram xjgamma.json \
  --histogram-receipt xjgamma.json.receipt.json \
  --dataset config/datasets/pp_engineering_fixture.yaml \
  --output xjgamma.png

photonjet plot emit-root-header \
  --histogram xjgamma.json \
  --histogram-receipt xjgamma.json.receipt.json \
  --dataset config/datasets/pp_engineering_fixture.yaml \
  --output PlotAnnotation.generated.h

photonjet plot verify-root-header \
  --histogram xjgamma.json \
  --histogram-receipt xjgamma.json.receipt.json \
  --dataset config/datasets/pp_engineering_fixture.yaml \
  --header PlotAnnotation.generated.h

photonjet skim \
  --input PhotonJetTrees_v1.root \
  --region A \
  --output recoil_pairs.root \
  --receipt recoil_pairs.json

photonjet purity \
  --input PhotonJetTrees_v1.root \
  --non-tight-definition bounded \
  --isolation-radius 0.4 \
  --output abcd_counts.json
```

The current public V1 starts at a validated `PhotonJetTrees_v1` file. The
mini-DST producer and direct-reference emitter remain an explicitly unproved
integration boundary; private migration adapters are not public prerequisites.

Each implemented stage is directly callable. The command only parses arguments
and dispatches to the documented Python API. The histogram receipt binds the
executed cuts and dataset manifest to the payload; the plotting commands
derive their labels from that receipt and expose no free-form cut-label
argument.
