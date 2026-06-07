#!/usr/bin/env python3
"""Build the OOF-vs-simple-holdout stack comparison report and slide."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


W, H = 2560, 1440
DPI = 200
COLORS = {
    "ink": "#151515",
    "muted": "#5f6368",
    "line": "#d7dbe0",
    "gray": "#eef1f4",
    "violet": "#6b4ea3",
    "teal": "#2f7f75",
    "green": "#3b8063",
    "amber": "#f4df96",
    "bluegray": "#526a7a",
    "white": "#ffffff",
}


def json_ready(value):
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, (np.floating, float)):
        f = float(value)
        return f if math.isfinite(f) else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_ready(payload), indent=2, sort_keys=True) + "\n")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size <= 0:
        raise SystemExit(f"Missing required compact artifact: {path}")


def load_domain(root: Path, domain: str, method: str) -> dict:
    droot = root / domain
    required = [
        "campaign_manifest.json",
        "partition_qa.json",
        "model_metrics.csv",
        "model_metrics.json",
        "stratified_metrics.csv",
        "leakage_qa.json",
        "score_correlations.json",
    ]
    for name in required:
        require_file(droot / name)
    manifest = read_json(droot / "campaign_manifest.json")
    leakage = read_json(droot / "leakage_qa.json")
    metrics = pd.read_csv(droot / "model_metrics.csv")
    strat = pd.read_csv(droot / "stratified_metrics.csv")
    return {
        "root": str(droot),
        "domain": domain,
        "method": method,
        "manifest": manifest,
        "leakage": leakage,
        "metrics": metrics,
        "stratified": strat,
    }


def locked_test_digest(payload: dict) -> str:
    return (
        payload["manifest"]
        .get("partition_qa", {})
        .get("partitions", {})
        .get("locked_test", {})
        .get("event_digest", {})
        .get("digest_blake2b16", "")
    )


def locked_test_rows(payload: dict) -> pd.DataFrame:
    df = payload["metrics"].copy()
    out = df[df["region"] == "locked_test"].copy()
    out["domain"] = payload["domain"]
    out["method"] = payload["method"]
    out["stack_training_mode"] = payload["manifest"].get("stack_training_mode", "oof5")
    return out


def metric_column(df: pd.DataFrame, suffix: str) -> str | None:
    for col in df.columns:
        if col.endswith(suffix):
            return col
    return None


def best_stack(df: pd.DataFrame) -> dict:
    stacks = df[~df["model"].isin(["BDT", "MLP"])].copy()
    if stacks.empty:
        fallback = df.sort_values("weighted_auc", ascending=False).iloc[0].to_dict()
        fallback["is_stack"] = False
        return fallback
    row = stacks.sort_values("weighted_auc", ascending=False).iloc[0].to_dict()
    row["is_stack"] = True
    return row


def model_value(df: pd.DataFrame, model: str, column: str) -> float:
    rows = df[df["model"] == model]
    if rows.empty or column not in rows:
        return math.nan
    return float(rows.iloc[0][column])


def label_model(name: str) -> str:
    labels = {
        "BDT": "BDT",
        "MLP": "MLP",
        "score_only_logistic": "Logistic stack",
        "score_only_gbm": "GBM stack",
        "score_only_mlp": "MLP stack",
        "score_context_logistic": "Logistic + context",
        "score_context_gbm": "GBM + context",
        "score_context_mlp": "MLP + context",
    }
    return labels.get(str(name), str(name).replace("_", " "))


def validate_pair(oof: dict, simple: dict) -> dict:
    domain = oof["domain"]
    oof_digest = locked_test_digest(oof)
    simple_digest = locked_test_digest(simple)
    oof_leak = oof["leakage"]
    simple_leak = simple["leakage"]
    return {
        "domain": domain,
        "locked_test_digest_oof": oof_digest,
        "locked_test_digest_simple_holdout": simple_digest,
        "locked_test_digest_match": bool(oof_digest and oof_digest == simple_digest),
        "oof_leakage_status": oof_leak.get("status"),
        "simple_holdout_leakage_status": simple_leak.get("status"),
        "oof_training_uses_source_sample_as_label": bool(
            oof["manifest"].get("class_definition_checks", {}).get("training_uses_source_sample_as_label", False)
        ),
        "simple_training_uses_source_sample_as_label": bool(
            simple["manifest"].get("class_definition_checks", {}).get("training_uses_source_sample_as_label", False)
        ),
        "simple_stack_mode": simple["manifest"].get("stack_training_mode"),
        "simple_stack_score_finite": bool(simple_leak.get("score_columns_finite", False)),
        "oof_stack_score_finite": bool(oof_leak.get("score_columns_finite", True)),
        "status": "pass",
    }


def build_summary_rows(payloads: dict[str, dict]) -> tuple[list[dict], list[dict]]:
    metric_frames = []
    summary_rows = []
    for domain in ("pp", "auau"):
        for method in ("current_oof", "simple_holdout"):
            metric_frames.append(locked_test_rows(payloads[f"{domain}:{method}"]))
    metrics = pd.concat(metric_frames, ignore_index=True)
    fake_col = metric_column(metrics, "weighted_background_fake_rate")
    for domain in ("pp", "auau"):
        domain_rows = metrics[metrics["domain"] == domain]
        for method in ("current_oof", "simple_holdout"):
            df = domain_rows[domain_rows["method"] == method]
            best = best_stack(df)
            bdt_auc = model_value(df, "BDT", "weighted_auc")
            mlp_auc = model_value(df, "MLP", "weighted_auc")
            best_auc = float(best.get("weighted_auc", math.nan))
            row = {
                "domain": domain,
                "method": method,
                "best_stack_model": best.get("model"),
                "best_stack_label": label_model(best.get("model")),
                "bdt_weighted_auc": bdt_auc,
                "mlp_weighted_auc": mlp_auc,
                "best_stack_weighted_auc": best_auc,
                "best_stack_minus_bdt_auc": best_auc - bdt_auc if math.isfinite(best_auc) and math.isfinite(bdt_auc) else math.nan,
                "best_stack_minus_mlp_auc": best_auc - mlp_auc if math.isfinite(best_auc) and math.isfinite(mlp_auc) else math.nan,
            }
            if fake_col:
                row["best_stack_wp80_weighted_fake_rate"] = float(best.get(fake_col, math.nan))
                row["bdt_wp80_weighted_fake_rate"] = model_value(df, "BDT", fake_col)
                row["mlp_wp80_weighted_fake_rate"] = model_value(df, "MLP", fake_col)
            summary_rows.append(row)
    comparison_rows = []
    for domain in ("pp", "auau"):
        oof = next(row for row in summary_rows if row["domain"] == domain and row["method"] == "current_oof")
        simple = next(row for row in summary_rows if row["domain"] == domain and row["method"] == "simple_holdout")
        delta = oof["best_stack_weighted_auc"] - simple["best_stack_weighted_auc"]
        favored = "neither"
        if math.isfinite(delta) and delta > 0.003:
            favored = "5-fold OOF"
        elif math.isfinite(delta) and delta < -0.003:
            favored = "simple holdout"
        comparison_rows.append(
            {
                "domain": domain,
                "current_oof_best_stack": oof["best_stack_model"],
                "simple_holdout_best_stack": simple["best_stack_model"],
                "current_oof_best_stack_weighted_auc": oof["best_stack_weighted_auc"],
                "simple_holdout_best_stack_weighted_auc": simple["best_stack_weighted_auc"],
                "oof_minus_simple_weighted_auc": delta,
                "favored_by_weighted_auc": favored,
            }
        )
    return metrics.to_dict("records"), summary_rows + comparison_rows


def method_takeaway(comparisons: list[dict]) -> str:
    favored = [row["favored_by_weighted_auc"] for row in comparisons if "favored_by_weighted_auc" in row]
    if favored and all(item == "5-fold OOF" for item in favored):
        return "Locked-test evidence favors 5-fold OOF in both domains."
    if favored and all(item == "simple holdout" for item in favored):
        return "Locked-test evidence favors simple holdout in both domains."
    if favored and "5-fold OOF" in favored and "simple holdout" not in favored:
        return "Locked-test evidence leans toward 5-fold OOF, but the gain is domain-dependent."
    if favored and "simple holdout" in favored and "5-fold OOF" not in favored:
        return "Locked-test evidence leans toward simple holdout, but the gain is domain-dependent."
    return "Locked-test evidence does not show a clear method-level advantage."


def setup_plot() -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": COLORS["line"],
            "axes.labelcolor": COLORS["ink"],
            "xtick.color": COLORS["ink"],
            "ytick.color": COLORS["ink"],
            "axes.titlecolor": COLORS["ink"],
        }
    )


def add_box(ax, x, y, w, h, fc, ec=None):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.008,rounding_size=0.01",
        transform=ax.transAxes,
        linewidth=1.1,
        edgecolor=ec or fc,
        facecolor=fc,
    )
    ax.add_patch(patch)
    return patch


def add_text(ax, x, y, text, size=18, weight="normal", color=None, ha="left", va="top", **kwargs):
    return ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        fontsize=size,
        fontweight=weight,
        color=color or COLORS["ink"],
        ha=ha,
        va=va,
        **kwargs,
    )


def plot_domain_bars(fig, rect, domain: str, rows: list[dict]) -> None:
    ax = fig.add_axes(rect)
    domain_rows = [row for row in rows if row.get("domain") == domain and "best_stack_weighted_auc" in row]
    oof = next(row for row in domain_rows if row["method"] == "current_oof")
    simple = next(row for row in domain_rows if row["method"] == "simple_holdout")
    labels = ["BDT", "MLP", "Simple\nbest stack", "OOF\nbest stack"]
    vals = [
        oof["bdt_weighted_auc"],
        oof["mlp_weighted_auc"],
        simple["best_stack_weighted_auc"],
        oof["best_stack_weighted_auc"],
    ]
    colors = [COLORS["bluegray"], COLORS["teal"], COLORS["violet"], COLORS["green"]]
    x = np.arange(len(vals))
    ax.bar(x, vals, color=colors, width=0.68)
    finite = [v for v in vals if math.isfinite(v)]
    if finite:
        ax.set_ylim(max(0.0, min(finite) - 0.035), min(1.0, max(finite) + 0.025))
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=15)
    ax.set_ylabel("Locked-test weighted AUC", fontsize=16)
    ax.set_title("pp" if domain == "pp" else "Au+Au", fontsize=23, fontweight="bold", pad=12)
    ax.grid(axis="y", color=COLORS["line"], linewidth=0.8)
    ax.tick_params(axis="y", labelsize=15)
    for i, val in enumerate(vals):
        label = f"{val:.3f}" if math.isfinite(val) else "NA"
        ax.text(i, val + 0.002 if math.isfinite(val) else 0.02, label, ha="center", va="bottom", fontsize=17)
    ax.text(
        0.5,
        -0.26,
        f"Simple: {label_model(simple['best_stack_model'])}   |   OOF: {label_model(oof['best_stack_model'])}",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=15,
        color=COLORS["muted"],
    )


def render_slide(outdir: Path, summary_rows: list[dict], validation_rows: list[dict]) -> Path:
    setup_plot()
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    add_text(ax, 0.045, 0.94, "Simple holdout vs 5-fold OOF stacking", size=34, weight="bold")
    add_text(
        ax,
        0.045,
        0.885,
        "Same base samples, features, labels, stackers, weights, and locked-test split.\nOnly the stack-training score source changes.",
        size=19,
        color=COLORS["muted"],
        linespacing=1.18,
    )

    add_box(ax, 0.055, 0.715, 0.275, 0.105, COLORS["gray"])
    add_text(ax, 0.075, 0.795, "Simple holdout", 22, "bold")
    add_text(ax, 0.075, 0.76, "One held-out trainval block\nscored by base models\nthat excluded it.", 15.8, linespacing=1.12)
    add_box(ax, 0.362, 0.715, 0.275, 0.105, "#ece7f5")
    add_text(ax, 0.382, 0.795, "5-fold OOF", 22, "bold")
    add_text(ax, 0.382, 0.76, "The same honesty rule\nrepeated over all\ntrainval folds.", 15.8, linespacing=1.12)
    add_box(ax, 0.67, 0.715, 0.275, 0.105, "#e6f1ef")
    add_text(ax, 0.69, 0.795, "Evaluation", 22, "bold")
    add_text(ax, 0.69, 0.76, "Both methods are judged\nonly on the same\nlocked test set.", 15.8, linespacing=1.12)

    plot_domain_bars(fig, [0.075, 0.335, 0.39, 0.255], "pp", summary_rows)
    plot_domain_bars(fig, [0.535, 0.335, 0.39, 0.255], "auau", summary_rows)

    comparisons = [row for row in summary_rows if "favored_by_weighted_auc" in row]
    digest_ok = all(row.get("locked_test_digest_match") for row in validation_rows)
    leak_ok = all(row.get("status") == "pass" for row in validation_rows)
    lines = []
    for row in comparisons:
        domain_label = "pp" if row["domain"] == "pp" else "Au+Au"
        lines.append(
            f"{domain_label}: OOF - simple = {row['oof_minus_simple_weighted_auc']:+.4f} weighted AUC "
            f"({row['favored_by_weighted_auc']})"
        )
    add_box(ax, 0.055, 0.045, 0.89, 0.16, COLORS["amber"], ec="#d2bd6c")
    add_text(ax, 0.075, 0.175, "Locked-test takeaway", 23, "bold")
    add_text(ax, 0.075, 0.135, method_takeaway(comparisons), 20, "bold")
    add_text(ax, 0.075, 0.10, "\n".join(lines), 16.5, color=COLORS["ink"], linespacing=1.16)
    add_text(
        ax,
        0.945,
        0.022,
        f"QA: locked-test digests {'match' if digest_ok else 'mismatch'}; leakage {'pass' if leak_ok else 'check'}",
        size=15,
        color=COLORS["muted"],
        ha="right",
    )

    path = outdir / "simple_holdout_vs_5fold_oof_stacking.png"
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def write_report(outdir: Path, summary_rows: list[dict], validation_rows: list[dict], args) -> Path:
    comparisons = [row for row in summary_rows if "favored_by_weighted_auc" in row]
    lines = [
        "# Simple holdout vs 5-fold OOF stacking",
        "",
        "## Contract",
        "",
        "- Simple holdout is the intuitive honest baseline: one trainval block is excluded from base BDT/MLP training, scored by those base models, then used to train the stackers.",
        "- 5-fold OOF is the data-efficient extension: every trainval fold receives base scores from base models that excluded that fold.",
        "- Both methods are evaluated on the same locked test split whenever the locked-test event digests match.",
        "- Training labels remain `is_signal == 1` vs `is_signal == 0`; `source_sample` is provenance and overlay/sample QA, not the supervised class label.",
        "",
        "## QA",
        "",
    ]
    for row in validation_rows:
        lines.append(
            f"- {row['domain']}: digest_match={row['locked_test_digest_match']}, "
            f"oof_leakage={row['oof_leakage_status']}, simple_leakage={row['simple_holdout_leakage_status']}, "
            f"source_as_label={row['oof_training_uses_source_sample_as_label'] or row['simple_training_uses_source_sample_as_label']}"
        )
    lines.extend(["", "## Locked-Test Comparison", ""])
    for row in comparisons:
        lines.append(
            f"- {row['domain']}: OOF best stack `{row['current_oof_best_stack']}` AUC "
            f"{row['current_oof_best_stack_weighted_auc']:.6f}; simple best stack "
            f"`{row['simple_holdout_best_stack']}` AUC {row['simple_holdout_best_stack_weighted_auc']:.6f}; "
            f"OOF-simple = {row['oof_minus_simple_weighted_auc']:+.6f}; favored = {row['favored_by_weighted_auc']}."
        )
    lines.extend(["", f"Bottom line: {method_takeaway(comparisons)}", ""])
    lines.extend(["## Inputs", "", f"- OOF root: `{args.oof_root}`", f"- Simple-holdout root: `{args.simple_holdout_root}`", ""])
    path = outdir / "comparison_report.md"
    path.write_text("\n".join(lines))
    return path


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--oof-root", type=Path, required=True, help="Local compact output root for current_oof with pp/ and auau/ subdirs.")
    ap.add_argument("--simple-holdout-root", type=Path, required=True, help="Local compact output root for simple_holdout with pp/ and auau/ subdirs.")
    ap.add_argument("--outdir", type=Path, required=True)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    payloads = {}
    for domain in ("pp", "auau"):
        payloads[f"{domain}:current_oof"] = load_domain(args.oof_root, domain, "current_oof")
        payloads[f"{domain}:simple_holdout"] = load_domain(args.simple_holdout_root, domain, "simple_holdout")

    validation_rows = [validate_pair(payloads[f"{domain}:current_oof"], payloads[f"{domain}:simple_holdout"]) for domain in ("pp", "auau")]
    for row in validation_rows:
        if not row["locked_test_digest_match"]:
            row["status"] = "fail"
        if row["oof_leakage_status"] != "pass" or row["simple_holdout_leakage_status"] != "pass":
            row["status"] = "fail"
        if row["oof_training_uses_source_sample_as_label"] or row["simple_training_uses_source_sample_as_label"]:
            row["status"] = "fail"
        if not row["simple_stack_score_finite"] or not row["oof_stack_score_finite"]:
            row["status"] = "fail"

    metrics_rows, summary_rows = build_summary_rows(payloads)
    write_csv(args.outdir / "locked_test_metrics.csv", metrics_rows)
    write_json(args.outdir / "locked_test_metrics.json", {"schema": "RJ_OOF_VS_SIMPLE_LOCKED_TEST_METRICS_V1", "rows": metrics_rows})
    write_csv(args.outdir / "method_summary.csv", summary_rows)
    write_json(args.outdir / "method_summary.json", {"schema": "RJ_OOF_VS_SIMPLE_METHOD_SUMMARY_V1", "rows": summary_rows})
    write_json(args.outdir / "leakage_validation.json", {"schema": "RJ_OOF_VS_SIMPLE_LEAKAGE_VALIDATION_V1", "rows": validation_rows})
    report = write_report(args.outdir, summary_rows, validation_rows, args)
    slide = render_slide(args.outdir, summary_rows, validation_rows)
    write_json(
        args.outdir / "comparison_manifest.json",
        {
            "schema": "RJ_OOF_VS_SIMPLE_COMPARISON_MANIFEST_V1",
            "oof_root": args.oof_root,
            "simple_holdout_root": args.simple_holdout_root,
            "outdir": args.outdir,
            "validation_status": "pass" if all(row["status"] == "pass" for row in validation_rows) else "fail",
            "artifacts": {
                "locked_test_metrics_csv": args.outdir / "locked_test_metrics.csv",
                "locked_test_metrics_json": args.outdir / "locked_test_metrics.json",
                "method_summary_csv": args.outdir / "method_summary.csv",
                "method_summary_json": args.outdir / "method_summary.json",
                "leakage_validation_json": args.outdir / "leakage_validation.json",
                "comparison_report_md": report,
                "slide_png": slide,
            },
        },
    )
    print(json.dumps({"status": "READY", "outdir": str(args.outdir), "slide": str(slide)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
