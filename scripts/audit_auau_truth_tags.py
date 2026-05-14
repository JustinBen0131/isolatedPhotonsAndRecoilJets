#!/usr/bin/env python3
"""Audit AuAu photon-ID truth signal/background labels in training trees.

The audit reads AuAuPhotonIDTrainingTree ROOT files directly so event identity
is available.  It checks whether embedded-photon samples mostly have one
signal-tagged reco photon per tagged event, and whether embedded-inclusive
samples have only a rare truth-signal tag rate.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REQUIRED_BRANCHES = ["run", "evt", "is_signal", "cluster_Et", "centrality"]
DEFAULT_PT_EDGES = [5.0, 10.0, 15.0, 20.0, 25.0, 35.0, 40.0]
DEFAULT_CENT_EDGES = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0]


@dataclass
class EventAgg:
    candidates: int = 0
    signal_tags: int = 0
    max_signal_et: float = math.nan
    min_signal_et: float = math.nan
    max_signal_cent: float = math.nan
    min_signal_cent: float = math.nan

    def add(self, is_signal: int, et: float, cent: float) -> None:
        self.candidates += 1
        if int(is_signal) == 1:
            self.signal_tags += 1
            if math.isfinite(et):
                if not math.isfinite(self.max_signal_et) or et > self.max_signal_et:
                    self.max_signal_et = et
                if not math.isfinite(self.min_signal_et) or et < self.min_signal_et:
                    self.min_signal_et = et
            if math.isfinite(cent):
                if not math.isfinite(self.max_signal_cent) or cent > self.max_signal_cent:
                    self.max_signal_cent = cent
                if not math.isfinite(self.min_signal_cent) or cent < self.min_signal_cent:
                    self.min_signal_cent = cent


def parse_edges(text: str, default: list[float]) -> list[float]:
    if not text:
        return list(default)
    vals = [float(x) for x in re.split(r"[,:;\s]+", text.strip()) if x]
    if len(vals) < 2 or any(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)):
        raise SystemExit(f"Invalid edge list: {text}")
    return vals


def bin_key(value: float, edges: list[float]) -> str | None:
    if not math.isfinite(value):
        return None
    for lo, hi in zip(edges[:-1], edges[1:]):
        if value >= lo and value < hi:
            return f"{lo:g}_{hi:g}"
    return None


def discover_paths(source: Path | None, manifest: Path | None) -> list[Path]:
    paths: list[Path] = []
    if manifest is not None:
        if not manifest.is_file():
            raise SystemExit(f"Manifest does not exist: {manifest}")
        paths.extend(Path(line.strip()) for line in manifest.read_text().splitlines()
                     if line.strip() and not line.strip().startswith("#"))
    elif source is not None:
        candidates = [
            source / "manifests" / "training_roots.list",
            source / "training_roots.list",
        ]
        for candidate in candidates:
            if candidate.is_file():
                return discover_paths(None, candidate)
        paths.extend(sorted(source.rglob("*.root")))
    else:
        raise SystemExit("Provide --source or --manifest")

    paths = [p for p in paths if str(p)]
    if not paths:
        raise SystemExit("No ROOT files found")
    return paths


def sample_from_path(path: Path) -> tuple[str, str]:
    text = str(path)
    lower = text.lower()
    if "embeddedphoton12" in lower:
        return "embeddedPhoton12", "signal_source"
    if "embeddedphoton20" in lower:
        return "embeddedPhoton20", "signal_source"
    if "embeddedjet12" in lower:
        return "embeddedJet12", "inclusive_source"
    if "embeddedjet20" in lower:
        return "embeddedJet20", "inclusive_source"
    if "/signal/" in lower or "signal:" in lower:
        return "signalUnknown", "signal_source"
    if "/background/" in lower or "inclusive" in lower or "jet" in lower:
        return "inclusiveUnknown", "inclusive_source"
    return "unknown", "unknown"


def limit_paths(paths: list[Path], max_files_per_sample: int, max_total_files: int) -> list[Path]:
    if max_files_per_sample <= 0 and max_total_files <= 0:
        return paths
    kept: list[Path] = []
    per_sample: Counter[str] = Counter()
    for path in paths:
        sample, _kind = sample_from_path(path)
        if max_files_per_sample > 0 and per_sample[sample] >= max_files_per_sample:
            continue
        kept.append(path)
        per_sample[sample] += 1
        if max_total_files > 0 and len(kept) >= max_total_files:
            break
    return kept


def init_bin_counts(edges: list[float]) -> dict[str, dict[str, int]]:
    return {
        f"{lo:g}_{hi:g}": {"candidates": 0, "signal_candidates": 0}
        for lo, hi in zip(edges[:-1], edges[1:])
    }


def safe_rate(num: int | float, den: int | float) -> float:
    return float(num) / float(den) if den else math.nan


def json_ready(obj):
    if isinstance(obj, dict):
        return {str(k): json_ready(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_ready(v) for v in obj]
    if isinstance(obj, tuple):
        return [json_ready(v) for v in obj]
    if isinstance(obj, (int, str, bool)) or obj is None:
        return obj
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    try:
        import numpy as np

        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            val = float(obj)
            return val if math.isfinite(val) else None
    except Exception:
        pass
    return str(obj)


def write_csv(path: Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def make_plots(outdir: Path, summary: dict, pt_rows: list[dict], cent_rows: list[dict]) -> list[str]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as exc:
        print(f"[truthTagAudit] plot generation skipped: {exc}", flush=True)
        return []

    plot_dir = outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    made: list[str] = []

    signal_mult = summary["by_source_class"].get("signal_source", {}).get("event_signal_multiplicity", {})
    if signal_mult:
        keys = sorted(int(k) for k in signal_mult)
        vals = [signal_mult.get(str(k), 0) for k in keys]
        fig, ax = plt.subplots(figsize=(6.8, 4.4))
        ax.bar(keys, vals, color="#0072B2")
        ax.set_xlabel("Signal-tagged reco photons per candidate-bearing event")
        ax.set_ylabel("Events")
        ax.set_title("Embedded-photon signal-tag multiplicity")
        ax.grid(axis="y", alpha=0.25)
        out = plot_dir / "signal_event_signal_multiplicity.png"
        fig.tight_layout()
        fig.savefig(out, dpi=180)
        plt.close(fig)
        made.append(str(out))

    def rate_plot(rows: list[dict], key_name: str, xlabel: str, out_name: str) -> None:
        incl = [r for r in rows if r.get("owner") == "inclusive_source"]
        if not incl:
            return
        labels = [str(r[key_name]).replace("_", "-") for r in incl]
        rates = [float(r["candidate_signal_fraction"]) for r in incl]
        counts = [int(r["candidates"]) for r in incl]
        x = np.arange(len(labels))
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        ax.bar(x, rates, color="#D55E00")
        ax.set_xticks(x, labels, rotation=35, ha="right")
        ax.set_ylabel("Signal-tagged candidate fraction")
        ax.set_xlabel(xlabel)
        ax.set_title("Embedded-inclusive truth-signal contamination")
        ax.grid(axis="y", alpha=0.25)
        for i, (rate, n) in enumerate(zip(rates, counts)):
            ax.text(i, rate, f"{rate:.3g}\nN={n}", ha="center", va="bottom", fontsize=8)
        out = plot_dir / out_name
        fig.tight_layout()
        fig.savefig(out, dpi=180)
        plt.close(fig)
        made.append(str(out))

    rate_plot(pt_rows, "pt_bin", r"cluster $E_T$ [GeV]", "inclusive_signal_tag_rate_by_et.png")
    rate_plot(cent_rows, "centrality_bin", "centrality [%]", "inclusive_signal_tag_rate_by_centrality.png")

    sample_rows = []
    for sample, payload in summary["by_sample"].items():
        sample_rows.append((sample, payload["signal_candidates"], payload["background_candidates"]))
    if sample_rows:
        sample_rows.sort()
        labels = [x[0] for x in sample_rows]
        sig = np.asarray([x[1] for x in sample_rows], dtype=float)
        bkg = np.asarray([x[2] for x in sample_rows], dtype=float)
        x = np.arange(len(labels))
        fig, ax = plt.subplots(figsize=(8.0, 4.6))
        ax.bar(x, bkg, color="#BBBBBB", label="is_signal = 0")
        ax.bar(x, sig, bottom=bkg, color="#0072B2", label="is_signal = 1")
        ax.set_yscale("log")
        ax.set_xticks(x, labels, rotation=25, ha="right")
        ax.set_ylabel("Candidates")
        ax.set_title("Candidate label composition by sample")
        ax.legend(frameon=False)
        ax.grid(axis="y", alpha=0.25)
        out = plot_dir / "candidate_label_composition_by_sample.png"
        fig.tight_layout()
        fig.savefig(out, dpi=180)
        plt.close(fig)
        made.append(str(out))

    return made


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=None)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    parser.add_argument("--max-files-per-sample", type=int, default=0)
    parser.add_argument("--max-total-files", type=int, default=0)
    parser.add_argument("--pt-edges", default="")
    parser.add_argument("--centrality-edges", default="")
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--inclusive-check-rate", type=float, default=0.03)
    parser.add_argument("--inclusive-fail-rate", type=float, default=0.10)
    parser.add_argument("--min-signal-single-tag-fraction", type=float, default=0.95)
    args = parser.parse_args()

    try:
        import numpy as np
        import uproot
    except Exception as exc:
        raise SystemExit("This audit needs uproot and numpy; use RJ_ML_PYTHON on SDCC.") from exc

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    pt_edges = parse_edges(args.pt_edges, DEFAULT_PT_EDGES)
    cent_edges = parse_edges(args.centrality_edges, DEFAULT_CENT_EDGES)

    paths = discover_paths(args.source, args.manifest)
    before_count = len(paths)
    paths = limit_paths(paths, args.max_files_per_sample, args.max_total_files)
    if not paths:
        raise SystemExit("All ROOT files were filtered out by limits")

    event_csv = outdir / "truth_tag_event_multiplicity.csv"
    dup_csv = outdir / "truth_tag_duplicate_events.csv"
    event_fields = [
        "sample", "source_class", "source_file", "run", "evt",
        "candidates", "signal_tags", "background_tags",
        "max_signal_et", "min_signal_et", "max_signal_cent", "min_signal_cent",
    ]
    duplicate_fields = event_fields

    by_sample: dict[str, dict] = {}
    by_class: dict[str, dict] = {}
    by_pt: dict[str, dict[str, dict[str, int]]] = defaultdict(lambda: init_bin_counts(pt_edges))
    by_cent: dict[str, dict[str, dict[str, int]]] = defaultdict(lambda: init_bin_counts(cent_edges))
    file_rows: list[dict] = []
    missing_rows: list[dict] = []
    duplicate_count = 0

    def ensure_bucket(container: dict[str, dict], key: str) -> dict:
        if key not in container:
            container[key] = {
                "files": 0,
                "events_with_candidates": 0,
                "events_with_signal": 0,
                "events_with_exactly_one_signal": 0,
                "events_with_multiple_signal": 0,
                "candidates": 0,
                "signal_candidates": 0,
                "background_candidates": 0,
                "event_signal_multiplicity": Counter(),
            }
        return container[key]

    with event_csv.open("w", newline="") as event_handle, dup_csv.open("w", newline="") as dup_handle:
        event_writer = csv.DictWriter(event_handle, fieldnames=event_fields)
        dup_writer = csv.DictWriter(dup_handle, fieldnames=duplicate_fields)
        event_writer.writeheader()
        dup_writer.writeheader()

        for idx, path in enumerate(paths, 1):
            sample, source_class = sample_from_path(path)
            if args.progress_every > 0 and (idx == 1 or idx % args.progress_every == 0 or idx == len(paths)):
                print(f"[truthTagAudit] {idx}/{len(paths)} {sample} {path}", flush=True)

            sample_bucket = ensure_bucket(by_sample, sample)
            class_bucket = ensure_bucket(by_class, source_class)
            sample_bucket["files"] += 1
            class_bucket["files"] += 1

            try:
                with uproot.open(path) as root_file:
                    if args.tree not in root_file:
                        missing_rows.append({"file": str(path), "reason": f"missing tree {args.tree}"})
                        continue
                    tree = root_file[args.tree]
                    keys = set(tree.keys())
                    missing = [name for name in REQUIRED_BRANCHES if name not in keys]
                    if missing:
                        missing_rows.append({"file": str(path), "reason": "missing branches: " + ",".join(missing)})
                        continue
                    arr = tree.arrays(REQUIRED_BRANCHES, library="np")
            except Exception as exc:
                missing_rows.append({"file": str(path), "reason": f"open/read failed: {exc}"})
                continue

            labels = np.asarray(arr["is_signal"], dtype=np.int32)
            valid = (labels == 0) | (labels == 1)
            labels = labels[valid]
            runs = np.asarray(arr["run"])[valid]
            evts = np.asarray(arr["evt"])[valid]
            ets = np.asarray(arr["cluster_Et"], dtype=np.float64)[valid]
            cents = np.asarray(arr["centrality"], dtype=np.float64)[valid]

            n_rows = int(len(labels))
            n_signal = int(np.sum(labels == 1))
            n_background = int(np.sum(labels == 0))
            for bucket in (sample_bucket, class_bucket):
                bucket["candidates"] += n_rows
                bucket["signal_candidates"] += n_signal
                bucket["background_candidates"] += n_background

            for label, et, cent in zip(labels, ets, cents):
                key = bin_key(float(et), pt_edges)
                if key is not None:
                    for owner in (sample, source_class):
                        by_pt[owner][key]["candidates"] += 1
                        by_pt[owner][key]["signal_candidates"] += int(label == 1)
                ckey = bin_key(float(cent), cent_edges)
                if ckey is not None:
                    for owner in (sample, source_class):
                        by_cent[owner][ckey]["candidates"] += 1
                        by_cent[owner][ckey]["signal_candidates"] += int(label == 1)

            events: dict[tuple[int, int], EventAgg] = defaultdict(EventAgg)
            for run, evt, label, et, cent in zip(runs, evts, labels, ets, cents):
                events[(int(run), int(evt))].add(int(label), float(et), float(cent))

            file_events_with_signal = 0
            file_events_multi = 0
            for (run, evt), agg in sorted(events.items()):
                row = {
                    "sample": sample,
                    "source_class": source_class,
                    "source_file": str(path),
                    "run": run,
                    "evt": evt,
                    "candidates": agg.candidates,
                    "signal_tags": agg.signal_tags,
                    "background_tags": agg.candidates - agg.signal_tags,
                    "max_signal_et": agg.max_signal_et,
                    "min_signal_et": agg.min_signal_et,
                    "max_signal_cent": agg.max_signal_cent,
                    "min_signal_cent": agg.min_signal_cent,
                }
                event_writer.writerow(row)
                for bucket in (sample_bucket, class_bucket):
                    bucket["events_with_candidates"] += 1
                    if agg.signal_tags > 0:
                        bucket["events_with_signal"] += 1
                    if agg.signal_tags == 1:
                        bucket["events_with_exactly_one_signal"] += 1
                    if agg.signal_tags > 1:
                        bucket["events_with_multiple_signal"] += 1
                    bucket["event_signal_multiplicity"][int(agg.signal_tags)] += 1
                if agg.signal_tags > 0:
                    file_events_with_signal += 1
                if agg.signal_tags > 1:
                    file_events_multi += 1
                    if source_class == "signal_source":
                        dup_writer.writerow(row)
                        duplicate_count += 1

            file_rows.append({
                "sample": sample,
                "source_class": source_class,
                "source_file": str(path),
                "candidates": n_rows,
                "signal_candidates": n_signal,
                "background_candidates": n_background,
                "events_with_candidates": len(events),
                "events_with_signal": file_events_with_signal,
                "events_with_multiple_signal": file_events_multi,
                "event_signal_tag_rate": safe_rate(file_events_with_signal, len(events)),
                "candidate_signal_fraction": safe_rate(n_signal, n_rows),
            })

    def finalize_bucket(payload: dict) -> dict:
        out = dict(payload)
        out["event_signal_multiplicity"] = {str(k): int(v) for k, v in sorted(payload["event_signal_multiplicity"].items())}
        out["event_signal_tag_rate"] = safe_rate(payload["events_with_signal"], payload["events_with_candidates"])
        out["candidate_signal_fraction"] = safe_rate(payload["signal_candidates"], payload["candidates"])
        out["single_signal_fraction_among_tagged_events"] = safe_rate(
            payload["events_with_exactly_one_signal"],
            payload["events_with_signal"],
        )
        out["multi_signal_fraction_among_candidate_events"] = safe_rate(
            payload["events_with_multiple_signal"],
            payload["events_with_candidates"],
        )
        return out

    sample_summary = {k: finalize_bucket(v) for k, v in sorted(by_sample.items())}
    class_summary = {k: finalize_bucket(v) for k, v in sorted(by_class.items())}

    def bin_rows(source: dict[str, dict[str, dict[str, int]]], axis_name: str) -> list[dict]:
        rows = []
        for owner, bins in sorted(source.items()):
            source_class = sample_from_path(Path(owner))[1]
            if owner in class_summary:
                source_class = owner
            elif owner in sample_summary:
                source_class = sample_from_path(Path(owner))[1]
            for key, counts in bins.items():
                rows.append({
                    "owner": owner,
                    "source_class": source_class if owner in class_summary else sample_from_path(Path(owner))[1],
                    axis_name: key,
                    "candidates": counts["candidates"],
                    "signal_candidates": counts["signal_candidates"],
                    "background_candidates": counts["candidates"] - counts["signal_candidates"],
                    "candidate_signal_fraction": safe_rate(counts["signal_candidates"], counts["candidates"]),
                })
        return rows

    pt_rows = bin_rows(by_pt, "pt_bin")
    cent_rows = bin_rows(by_cent, "centrality_bin")

    write_csv(outdir / "truth_tag_file_summary.csv", file_rows, [
        "sample", "source_class", "source_file", "candidates", "signal_candidates",
        "background_candidates", "events_with_candidates", "events_with_signal",
        "events_with_multiple_signal", "event_signal_tag_rate", "candidate_signal_fraction",
    ])
    write_csv(outdir / "truth_tag_pt_bins.csv", pt_rows, [
        "owner", "source_class", "pt_bin", "candidates", "signal_candidates",
        "background_candidates", "candidate_signal_fraction",
    ])
    write_csv(outdir / "truth_tag_centrality_bins.csv", cent_rows, [
        "owner", "source_class", "centrality_bin", "candidates", "signal_candidates",
        "background_candidates", "candidate_signal_fraction",
    ])
    if missing_rows:
        write_csv(outdir / "truth_tag_missing_or_bad_files.csv", missing_rows, ["file", "reason"])

    notes: list[str] = [
        "Event denominator is candidate-bearing events in AuAuPhotonIDTrainingTree, not all generated/processed events.",
        "Event key is (sample, source_file, run, evt), so chunk-local event counters do not collide.",
    ]
    status = "PASS"

    signal_payload = class_summary.get("signal_source", {})
    inclusive_payload = class_summary.get("inclusive_source", {})
    signal_single = signal_payload.get("single_signal_fraction_among_tagged_events", math.nan)
    inclusive_rate = inclusive_payload.get("event_signal_tag_rate", math.nan)

    if missing_rows:
        status = "CHECK"
        notes.append(f"{len(missing_rows)} files were missing/read-bad; see truth_tag_missing_or_bad_files.csv.")
    if not signal_payload or signal_payload.get("signal_candidates", 0) <= 0:
        status = "FAIL"
        notes.append("No signal-tagged candidates found in embedded-photon source files.")
    elif math.isfinite(signal_single) and signal_single < args.min_signal_single_tag_fraction:
        status = "CHECK"
        notes.append(
            f"Signal-source single-tag fraction among tagged events is {signal_single:.4g}, "
            f"below {args.min_signal_single_tag_fraction:.4g}."
        )
    if not inclusive_payload:
        status = "CHECK" if status == "PASS" else status
        notes.append("No embedded-inclusive source files found.")
    elif math.isfinite(inclusive_rate):
        if inclusive_rate >= args.inclusive_fail_rate:
            status = "FAIL"
            notes.append(f"Inclusive-source event signal-tag rate is order-one-like: {inclusive_rate:.4g}.")
        elif inclusive_rate >= args.inclusive_check_rate:
            status = "CHECK" if status == "PASS" else status
            notes.append(
                f"Inclusive-source event signal-tag rate is {inclusive_rate:.4g}, "
                f"above check threshold {args.inclusive_check_rate:.4g}."
            )

    summary = {
        "schema": "RJ_AUAU_TRUTH_TAG_AUDIT_V1",
        "status": status,
        "source": str(args.source) if args.source else None,
        "manifest": str(args.manifest) if args.manifest else None,
        "tree": args.tree,
        "input_files_before_limits": before_count,
        "input_files_audited": len(paths),
        "max_files_per_sample": args.max_files_per_sample,
        "max_total_files": args.max_total_files,
        "pt_edges": pt_edges,
        "centrality_edges": cent_edges,
        "denominator": "candidate-bearing events in AuAuPhotonIDTrainingTree",
        "duplicate_signal_source_events": duplicate_count,
        "thresholds": {
            "inclusive_check_rate": args.inclusive_check_rate,
            "inclusive_fail_rate": args.inclusive_fail_rate,
            "min_signal_single_tag_fraction": args.min_signal_single_tag_fraction,
        },
        "by_sample": sample_summary,
        "by_source_class": class_summary,
        "bad_files": missing_rows,
        "notes": notes,
        "outputs": {
            "event_multiplicity_csv": str(event_csv),
            "duplicate_events_csv": str(dup_csv),
            "file_summary_csv": str(outdir / "truth_tag_file_summary.csv"),
            "pt_bins_csv": str(outdir / "truth_tag_pt_bins.csv"),
            "centrality_bins_csv": str(outdir / "truth_tag_centrality_bins.csv"),
        },
    }

    if not args.no_plots:
        plots = make_plots(outdir, summary, pt_rows, cent_rows)
        summary["outputs"]["plots"] = plots

    (outdir / "truth_tag_audit_summary.json").write_text(
        json.dumps(json_ready(summary), indent=2, sort_keys=True) + "\n"
    )
    print(
        "[truthTagAudit] "
        f"status={status} files={len(paths)} "
        f"signal_single_tag_fraction={signal_single if math.isfinite(signal_single) else 'nan'} "
        f"inclusive_event_signal_rate={inclusive_rate if math.isfinite(inclusive_rate) else 'nan'} "
        f"outdir={outdir}",
        flush=True,
    )
    return 0 if status in {"PASS", "CHECK"} else 3


if __name__ == "__main__":
    raise SystemExit(main())
