#!/usr/bin/env python3
"""Aggregate XGBoost split-gain usage for a routed AuAu BDT product."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


PRODUCT_RE = re.compile(
    r"auau_tight_bdt_(?P<product>.+)_pt_(?P<ptlo>\d{3})_(?P<pthi>\d{3})_cent_(?P<clo>\d{3})_(?P<chi>\d{3})_tmva\.metadata\.json$"
)

FEATURE_LABELS = {
    "cluster_Et": r"cluster $E_T$",
    "cluster_weta_cogx": r"$w_{\eta}$",
    "cluster_wphi_cogx": r"$w_{\phi}$",
    "cluster_weta33_cogx": r"$w_{\eta}^{3x3}$",
    "cluster_wphi33_cogx": r"$w_{\phi}^{3x3}$",
    "vertexz": r"$z_{\mathrm{vtx}}$",
    "cluster_Eta": r"cluster $\eta$",
    "e11_over_e33": r"$E_{11}/E_{33}$",
    "cluster_et1": "cluster et1",
    "cluster_et2": "cluster et2",
    "cluster_et3": "cluster et3",
    "cluster_et4": "cluster et4",
    "e32_over_e35": r"$E_{32}/E_{35}$",
    "cluster_weta35_cogx": r"$w_{\eta}^{3x5}$",
    "cluster_wphi53_cogx": r"$w_{\phi}^{5x3}$",
    "cluster_w32": r"$w_{32}$",
    "cluster_w52": r"$w_{52}$",
    "cluster_w72": r"$w_{72}$",
    "e11_over_e22": r"$E_{11}/E_{22}$",
    "e11_over_e13": r"$E_{11}/E_{13}$",
    "e11_over_e15": r"$E_{11}/E_{15}$",
    "e11_over_e17": r"$E_{11}/E_{17}$",
    "e11_over_e31": r"$E_{11}/E_{31}$",
    "e11_over_e51": r"$E_{11}/E_{51}$",
    "e11_over_e71": r"$E_{11}/E_{71}$",
    "e22_over_e33": r"$E_{22}/E_{33}$",
    "e22_over_e35": r"$E_{22}/E_{35}$",
    "e22_over_e37": r"$E_{22}/E_{37}$",
    "e22_over_e53": r"$E_{22}/E_{53}$",
    "cluster_weta_over_wphi": r"$w_{\eta}/w_{\phi}$",
    "cluster_weta33_over_wphi33": r"$w_{\eta}^{3x3}/w_{\phi}^{3x3}$",
    "centrality": "centrality",
    "reco_eiso_clip30": r"$E_T^{iso}$ clipped",
    "reco_eiso_over_cluster_Et": r"$E_T^{iso}/E_T^{cluster}$",
    "reco_eiso_signed_log1p": r"signed log$(1+E_T^{iso})$",
    "reco_eiso_r30": r"R=0.3 $E_T^{iso}$",
    "reco_eiso_r40": r"R=0.4 $E_T^{iso}$",
}

FEATURE_FAMILIES = {
    "kinematic/context": {
        "cluster_Et",
        "cluster_Eta",
        "vertexz",
        "centrality",
    },
    "shower widths": {
        "cluster_weta_cogx",
        "cluster_wphi_cogx",
        "cluster_weta33_cogx",
        "cluster_wphi33_cogx",
        "cluster_weta35_cogx",
        "cluster_wphi53_cogx",
        "cluster_w32",
        "cluster_w52",
        "cluster_w72",
        "cluster_weta_over_wphi",
        "cluster_weta33_over_wphi33",
    },
    "energy sharing": {
        "cluster_et1",
        "cluster_et2",
        "cluster_et3",
        "cluster_et4",
        "e11_over_e33",
        "e32_over_e35",
        "e11_over_e22",
        "e11_over_e13",
        "e11_over_e15",
        "e11_over_e17",
        "e11_over_e31",
        "e11_over_e51",
        "e11_over_e71",
        "e22_over_e33",
        "e22_over_e35",
        "e22_over_e37",
        "e22_over_e53",
    },
    "isolation inputs": {
        "reco_eiso_clip30",
        "reco_eiso_over_cluster_Et",
        "reco_eiso_signed_log1p",
        "reco_eiso_r30",
        "reco_eiso_r40",
    },
}

FAMILY_COLORS = {
    "kinematic/context": "#6B7280",
    "shower widths": "#0072B2",
    "energy sharing": "#009E73",
    "isolation inputs": "#CC79A7",
    "other": "#999999",
}


def feature_family(feature: str) -> str:
    for family, features in FEATURE_FAMILIES.items():
        if feature in features:
            return family
    return "other"


def parse_route(path: Path) -> tuple[str, float, float, float, float]:
    match = PRODUCT_RE.search(path.name)
    if not match:
        raise ValueError(f"Could not parse product/route from {path.name}")
    product = match.group("product")
    values = tuple(float(match.group(k)) for k in ("ptlo", "pthi", "clo", "chi"))
    return product, values[0], values[1], values[2], values[3]


def gain_and_counts_by_feature(xgb_json: Path, features: list[str]) -> tuple[dict[str, float], dict[str, int]]:
    data = json.loads(xgb_json.read_text())
    trees = data["learner"]["gradient_booster"]["model"]["trees"]
    gains: dict[str, float] = defaultdict(float)
    counts: dict[str, int] = defaultdict(int)
    for tree in trees:
        for left, right, idx, loss in zip(
            tree["left_children"],
            tree["right_children"],
            tree["split_indices"],
            tree["loss_changes"],
        ):
            if left < 0 and right < 0:
                continue
            if int(idx) < len(features):
                feature = features[int(idx)]
                gains[feature] += max(float(loss), 0.0)
                counts[feature] += 1
    return dict(gains), dict(counts)


def collect(model_dir: Path, product: str) -> tuple[list[dict], list[str], dict[str, object]]:
    rows: list[dict] = []
    feature_order: list[str] | None = None
    route_count = 0
    auc_values: list[float] = []
    for meta_path in sorted(model_dir.glob(f"auau_tight_bdt_{product}_pt_*_cent_*_tmva.metadata.json")):
        xgb_path = meta_path.with_name(meta_path.name.replace(".metadata.json", ".xgb.json"))
        if not xgb_path.exists():
            continue
        parsed_product, ptlo, pthi, clo, chi = parse_route(meta_path)
        if parsed_product != product:
            continue
        meta = json.loads(meta_path.read_text())
        features = list(meta["features"])
        if feature_order is None:
            feature_order = features
        if feature_order != features:
            raise SystemExit(f"Feature-order mismatch in {meta_path}")
        gains, counts = gain_and_counts_by_feature(xgb_path, features)
        total_gain = sum(gains.values())
        route_count += 1
        try:
            auc_values.append(float(meta["auc"]))
        except Exception:
            pass
        for feature in features:
            gain = gains.get(feature, 0.0)
            rows.append(
                {
                    "product": product,
                    "route": f"{ptlo:g}-{pthi:g} GeV, {clo:g}-{chi:g}%",
                    "pt_lo": ptlo,
                    "pt_hi": pthi,
                    "cent_lo": clo,
                    "cent_hi": chi,
                    "feature": feature,
                    "feature_label": FEATURE_LABELS.get(feature, feature),
                    "feature_family": feature_family(feature),
                    "route_gain": gain,
                    "route_gain_fraction": gain / total_gain if total_gain > 0 else math.nan,
                    "route_split_count": counts.get(feature, 0),
                    "route_auc": float(meta.get("auc", "nan")),
                }
            )
    if not rows:
        raise SystemExit(f"No metadata/xgb pairs found for product={product} in {model_dir}")
    meta_summary: dict[str, object] = {
        "n_routes": route_count,
        "n_features": len(feature_order or []),
        "route_auc_mean": float(sum(auc_values) / len(auc_values)) if auc_values else math.nan,
        "route_auc_min": float(min(auc_values)) if auc_values else math.nan,
        "route_auc_max": float(max(auc_values)) if auc_values else math.nan,
    }
    return rows, feature_order or [], meta_summary


def aggregate_rows(rows: list[dict], features: list[str]) -> list[dict]:
    total_gain = sum(float(r["route_gain"]) for r in rows)
    total_splits = sum(int(r["route_split_count"]) for r in rows)
    out: list[dict] = []
    for feature in features:
        gain = sum(float(r["route_gain"]) for r in rows if r["feature"] == feature)
        count = sum(int(r["route_split_count"]) for r in rows if r["feature"] == feature)
        out.append(
            {
                "feature": feature,
                "feature_label": FEATURE_LABELS.get(feature, feature),
                "feature_family": feature_family(feature),
                "split_gain": gain,
                "split_gain_fraction": gain / total_gain if total_gain > 0 else math.nan,
                "split_count": count,
                "split_count_fraction": count / total_splits if total_splits > 0 else math.nan,
            }
        )
    out.sort(key=lambda r: (float(r["split_gain_fraction"]), int(r["split_count"])), reverse=True)
    return out


def draw_plot(rows: list[dict], summary: dict[str, object], out: Path) -> None:
    values = list(rows)
    n = len(values)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig_h = max(12.0, 0.36 * n + 3.4)
    fig, ax = plt.subplots(figsize=(15.5, fig_h), dpi=185)
    fig.patch.set_facecolor("white")
    y = list(range(n))
    frac = [100.0 * float(r["split_gain_fraction"]) for r in values]
    labels = [str(r["feature_label"]) for r in values]
    colors = [FAMILY_COLORS.get(str(r["feature_family"]), "#999999") for r in values]

    ax.barh(y, frac, color=colors, height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=12.4)
    ax.invert_yaxis()
    ax.set_xlabel("Fraction of total split gain [%]", fontsize=15)
    ax.tick_params(axis="x", labelsize=12)
    ax.grid(axis="x", color="#D1D5DB", linewidth=1.0, alpha=0.7)
    ax.set_axisbelow(True)
    xmax = max(frac) * 1.18 if frac else 1.0
    ax.set_xlim(0.0, xmax)
    for yy, val, item in zip(y, frac, values):
        ax.text(
            val + xmax * 0.007,
            yy,
            f"{val:.1f}% ({int(item['split_count'])})",
            va="center",
            ha="left",
            fontsize=11.5,
            fontweight="bold" if val >= 2.0 else "normal",
        )

    legend_items = [
        Line2D([0], [0], color=FAMILY_COLORS[k], lw=8, label=k)
        for k in ("kinematic/context", "shower widths", "energy sharing", "isolation inputs")
    ]
    ax.legend(
        handles=legend_items,
        loc="lower right",
        bbox_to_anchor=(0.995, 0.015),
        frameon=True,
        facecolor="white",
        edgecolor="#E5E7EB",
        fontsize=11.8,
        title="Feature family",
        title_fontsize=12.4,
    )

    fig.text(0.067, 0.973, "Split-gain usage: binned BDT with isolation inputs", ha="left", va="top", fontsize=24, fontweight="bold")
    fig.text(
        0.067,
        0.945,
        f"8 $E_T$ x 7 centrality routed XGBoost BDTs; gain summed over {summary['n_routes']} route-specific models; labels show gain fraction and split count",
        ha="left",
        va="top",
        fontsize=13.3,
        color="#374151",
    )
    fig.text(0.067, 0.915, "sPHENIX", ha="left", va="top", fontsize=18, fontstyle="italic", fontweight="bold")
    fig.text(0.143, 0.915, " Internal", ha="left", va="top", fontsize=18)
    fig.text(
        0.067,
        0.890,
        "Photon12+20 signal vs Jet12+20+30 background, 15 < cluster $E_T$ < 35 GeV",
        ha="left",
        va="top",
        fontsize=12.0,
        color="#111827",
    )
    fig.tight_layout(rect=[0.05, 0.045, 0.985, 0.865])
    fig.savefig(out)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", type=Path, required=True)
    ap.add_argument("--product", default="globalEtCent1535_bdt_iso_ptCent7")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--tag", default="binned_bdt_iso_ptcent7_split_gain_all_features")
    args = ap.parse_args()

    route_rows, features, summary = collect(args.model_dir, args.product)
    agg_rows = aggregate_rows(route_rows, features)
    args.outdir.mkdir(parents=True, exist_ok=True)
    route_csv = args.outdir / f"{args.tag}_by_route.csv"
    agg_csv = args.outdir / f"{args.tag}.csv"
    png = args.outdir / f"{args.tag}.png"
    manifest = args.outdir / f"{args.tag}.json"

    with route_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(route_rows[0]))
        writer.writeheader()
        writer.writerows(route_rows)
    with agg_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(agg_rows[0]))
        writer.writeheader()
        writer.writerows(agg_rows)
    draw_plot(agg_rows, summary, png)
    manifest.write_text(
        json.dumps(
            {
                "schema": "AUAU_BINNED_BDT_SPLIT_GAIN_BREAKDOWN_V1",
                "model_dir": str(args.model_dir),
                "product": args.product,
                **summary,
                "outputs": {
                    "png": str(png),
                    "aggregate_csv": str(agg_csv),
                    "by_route_csv": str(route_csv),
                    "json": str(manifest),
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(png)
    print(agg_csv)
    print(route_csv)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
