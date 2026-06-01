# SDCC Scripts

Canonical source for scripts that may run on SDCC. Runtime anchors live under `runtime/`, campaign drivers under `workflows/`, pipeline entrypoints under `pipelines/`, transfer helpers under `transfer/`, and local ROOT merge helpers under `local_merge/root/`.

See `TRANSFER_MAP.tsv`, `SDCC_RUNTIME_INDEX.yaml`, and `IO_RUNTIME_INDEX.yaml` before uploading, moving, or changing output paths. Use `runtime/audit/checkout_hygiene_probe.py` for read-only checkout hygiene checks and `runtime/io/resolve_io_contract.py` to resolve canonical vs legacy campaign I/O paths.
