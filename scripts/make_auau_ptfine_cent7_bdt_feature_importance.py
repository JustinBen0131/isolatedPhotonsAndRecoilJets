#!/usr/bin/env python3
"""Plot XGBoost split-gain usage for the routed AuAu pT x cent7 BDT."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


FEATURE_LABELS = {
    "cluster_Et": "cluster_Et",
    "cluster_weta_cogx": "cluster_weta_cogx",
    "cluster_wphi_cogx": "cluster_wphi_cogx",
    "cluster_weta33_cogx": "cluster_weta33_cogx",
    "cluster_wphi33_cogx": "cluster_wphi33_cogx",
    "vertexz": "vertexz",
    "cluster_Eta": "cluster_Eta",
    "e11_over_e33": "e11_over_e33",
    "cluster_et1": "cluster_et1",
    "cluster_et2": "cluster_et2",
    "cluster_et3": "cluster_et3",
    "cluster_et4": "cluster_et4",
    "e32_over_e35": "e32_over_e35",
}

COLORS = {
    "cluster_Et": "#4C78A8",
    "cluster_weta_cogx": "#1F33CC",
    "cluster_wphi_cogx": "#36A2EB",
    "cluster_weta33_cogx": "#FF7F0E",
    "cluster_wphi33_cogx": "#D62728",
    "vertexz": "#8C564B",
    "cluster_Eta": "#9467BD",
    "e11_over_e33": "#2CA02C",
    "cluster_et1": "#17BECF",
    "cluster_et2": "#BCBD22",
    "cluster_et3": "#E377C2",
    "cluster_et4": "#7F7F7F",
    "e32_over_e35": "#111111",
}

ROUTE_RE = re.compile(
    r"ptFine_cent7_pt_(?P<ptlo>\d{3})_(?P<pthi>\d{3})_cent_(?P<clo>\d{3})_(?P<chi>\d{3})"
)


def font(name: str, size: int) -> ImageFont.FreeTypeFont:
    base = Path("/System/Library/Fonts/Supplemental")
    paths = {
        "regular": base / "Times New Roman.ttf",
        "bold": base / "Times New Roman Bold.ttf",
        "italic": base / "Times New Roman Italic.ttf",
        "bolditalic": base / "Times New Roman Bold Italic.ttf",
    }
    try:
        return ImageFont.truetype(str(paths[name]), size=size)
    except Exception:
        return ImageFont.load_default()


def parse_route(path: Path) -> tuple[float, float, float, float]:
    match = ROUTE_RE.search(path.stem)
    if not match:
        raise ValueError(f"Could not parse route from {path.name}")
    return tuple(float(match.group(k)) for k in ("ptlo", "pthi", "clo", "chi"))  # type: ignore[return-value]


def gain_by_feature(xgb_json: Path, features: list[str]) -> dict[str, float]:
    data = json.loads(xgb_json.read_text())
    trees = data["learner"]["gradient_booster"]["model"]["trees"]
    gain = defaultdict(float)
    for tree in trees:
        for left, right, idx, loss in zip(
            tree["left_children"],
            tree["right_children"],
            tree["split_indices"],
            tree["loss_changes"],
        ):
            if left < 0 and right < 0:
                continue
            if idx < len(features):
                gain[features[idx]] += max(float(loss), 0.0)
    return dict(gain)


def collect(model_dir: Path) -> tuple[list[dict], list[str]]:
    rows = []
    feature_order: list[str] | None = None
    for meta_path in sorted(model_dir.glob("auau_tight_bdt_ptFine_cent7_*_tmva.metadata.json")):
        xgb_path = meta_path.with_name(meta_path.name.replace(".metadata.json", ".xgb.json"))
        if not xgb_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        features = list(meta["features"])
        if feature_order is None:
            feature_order = features
        ptlo, pthi, clo, chi = parse_route(meta_path)
        gains = gain_by_feature(xgb_path, features)
        total = sum(gains.values())
        for feature in features:
            val = gains.get(feature, 0.0)
            rows.append(
                {
                    "route": f"{ptlo:g}-{pthi:g} GeV, {clo:g}-{chi:g}%",
                    "pt_lo": ptlo,
                    "pt_hi": pthi,
                    "cent_lo": clo,
                    "cent_hi": chi,
                    "feature": feature,
                    "feature_label": FEATURE_LABELS.get(feature, feature),
                    "gain": val,
                    "gain_fraction": val / total if total > 0 else math.nan,
                    "auc": float(meta.get("auc", "nan")),
                }
            )
    return rows, feature_order or []


def draw_sphenix(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    draw.text((x, y), "sPHENIX", fill=(0, 0, 0), font=font("bolditalic", 22))
    draw.text((x + 104, y + 1), "Internal", fill=(0, 0, 0), font=font("regular", 21))
    draw.text((x, y + 28), "Au+Au embedded photon-ID training", fill=(0, 0, 0), font=font("regular", 18))


def _text_width(draw: ImageDraw.ImageDraw, text: str, text_font: ImageFont.FreeTypeFont) -> float:
    bbox = draw.textbbox((0, 0), text, font=text_font)
    return bbox[2] - bbox[0]


def draw_et_text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    prefix: str,
    suffix: str,
    *,
    main_font: ImageFont.FreeTypeFont,
    sub_font: ImageFont.FreeTypeFont,
    fill: tuple[int, int, int] = (0, 0, 0),
    anchor: str = "mm",
) -> None:
    """Draw centered text containing E with a subscript T."""
    x, y = xy
    e_w = _text_width(draw, "E", main_font)
    t_w = _text_width(draw, "T", sub_font)
    pieces = [
        (prefix, main_font, 0.0),
        ("E", main_font, 0.0),
        ("T", sub_font, 0.30),
        (suffix, main_font, 0.0),
    ]
    total = _text_width(draw, prefix, main_font) + e_w + t_w + _text_width(draw, suffix, main_font)
    start_x = x - total / 2 if anchor == "mm" else x
    cursor = start_x
    for text, text_font, y_shift in pieces:
        draw.text((cursor, y + y_shift * text_font.size), text, fill=fill, font=text_font, anchor="lm")
        cursor += _text_width(draw, text, text_font)


def draw_overall(rows: list[dict], features: list[str], out: Path) -> None:
    total_gain = sum(r["gain"] for r in rows)
    values = [(f, sum(r["gain"] for r in rows if r["feature"] == f) / total_gain) for f in features]
    values = sorted(values, key=lambda x: x[1], reverse=True)
    w, h = 1700, 900
    img = Image.new("RGB", (w, h), "white")
    draw = ImageDraw.Draw(img)
    draw_et_text(
        draw,
        (w // 2, 34),
        "Split-gain usage across the routed ",
        " x cent7 BDT",
        main_font=font("bold", 36),
        sub_font=font("bold", 24),
    )
    draw.text((w // 2, 74), "56 bin-specific XGBoost models; gain summed over all routes and normalized", anchor="mm", fill=(85, 85, 85), font=font("regular", 22))
    draw_sphenix(draw, 84, 36)

    left, top, right, row_h = 390, 145, 1320, 48
    max_v = max(v for _, v in values) * 1.08
    draw.line((left, top + len(values) * row_h + 12, right, top + len(values) * row_h + 12), fill=(80, 80, 80), width=2)
    for i in range(6):
        x = left + (right - left) * i / 5
        val = max_v * i / 5
        draw.line((x, top - 6, x, top + len(values) * row_h + 12), fill=(230, 230, 230), width=1)
        draw.text((x, top + len(values) * row_h + 24), f"{100*val:.0f}", anchor="ma", fill=(50, 50, 50), font=font("regular", 18))
    draw.text(((left + right) // 2, h - 42), "Fraction of total split gain [%]", anchor="mm", fill=(0, 0, 0), font=font("regular", 24))

    for i, (feature, val) in enumerate(values):
        y = top + i * row_h
        label = FEATURE_LABELS.get(feature, feature)
        draw.text((left - 18, y + row_h / 2), label, anchor="rm", fill=(0, 0, 0), font=font("regular", 20))
        x1 = left + int((right - left) * val / max_v)
        draw.rounded_rectangle((left, y + 8, x1, y + row_h - 8), radius=5, fill=COLORS.get(feature, "#888888"))
        draw.text((x1 + 10, y + row_h / 2), f"{100*val:.1f}", anchor="lm", fill=(0, 0, 0), font=font("bold", 20))

    out.parent.mkdir(parents=True, exist_ok=True)
    img.save(out)


def draw_stacked_by_pt(rows: list[dict], features: list[str], out: Path) -> None:
    pt_bins = sorted({(r["pt_lo"], r["pt_hi"]) for r in rows})
    by_pt_feature = defaultdict(float)
    by_pt_total = defaultdict(float)
    for r in rows:
        key = (r["pt_lo"], r["pt_hi"])
        by_pt_feature[(key, r["feature"])] += r["gain"]
        by_pt_total[key] += r["gain"]
    w, h = 1500, 900
    img = Image.new("RGB", (w, h), "white")
    draw = ImageDraw.Draw(img)
    draw_et_text(
        draw,
        (w // 2, 34),
        "Feature mix used by the ",
        " x cent7 routed BDT",
        main_font=font("bold", 36),
        sub_font=font("bold", 24),
    )
    draw_et_text(
        draw,
        (w // 2, 76),
        "Each bar sums the seven centrality-specific BDTs in one ",
        " bin",
        main_font=font("regular", 23),
        sub_font=font("regular", 16),
        fill=(85, 85, 85),
    )
    draw_sphenix(draw, 84, 36)

    left, top, bottom = 135, 155, 710
    bar_w, gap = 105, 24
    plot_h = bottom - top
    y_max = 1.0
    for i in range(6):
        y = bottom - plot_h * i / 5
        val = y_max * i / 5
        draw.line((left - 8, y, left + len(pt_bins) * (bar_w + gap), y), fill=(230, 230, 230), width=1)
        draw.text((left - 18, y), f"{100*val:.0f}", anchor="rm", fill=(50, 50, 50), font=font("regular", 20))
    ylab = Image.new("RGBA", (360, 42), (255, 255, 255, 0))
    ydraw = ImageDraw.Draw(ylab)
    ydraw.text((180, 21), "Fraction of split gain [%]", anchor="mm", fill=(0, 0, 0), font=font("regular", 25))
    img.paste(ylab.rotate(90, expand=True), (12, 305), ylab.rotate(90, expand=True))

    for ip, pt in enumerate(pt_bins):
        x0 = left + ip * (bar_w + gap)
        y0 = bottom
        total = by_pt_total[pt]
        ordered = sorted(features, key=lambda f: by_pt_feature[(pt, f)], reverse=True)
        for feature in ordered:
            frac = by_pt_feature[(pt, feature)] / total if total > 0 else 0
            if frac <= 0:
                continue
            y1 = y0 - plot_h * frac / y_max
            draw.rectangle((x0, y1, x0 + bar_w, y0), fill=COLORS.get(feature, "#888888"), outline="white")
            y0 = y1
        draw.rectangle((x0, top, x0 + bar_w, bottom), outline=(80, 80, 80), width=1)
        draw.text((x0 + bar_w / 2, bottom + 18), f"{pt[0]:.0f}-{pt[1]:.0f}", anchor="ma", fill=(0, 0, 0), font=font("regular", 20))
    draw_et_text(
        draw,
        (left + len(pt_bins) * (bar_w + gap) / 2 - gap / 2, h - 58),
        "Candidate ",
        " bin [GeV]",
        main_font=font("regular", 27),
        sub_font=font("regular", 18),
    )

    leg_x, leg_y = 1215, 155
    draw.text((leg_x, leg_y - 40), "Feature", fill=(0, 0, 0), font=font("bold", 25))
    for i, feature in enumerate(features):
        y = leg_y + i * 34
        draw.rectangle((leg_x, y, leg_x + 23, y + 23), fill=COLORS.get(feature, "#888888"))
        draw.text((leg_x + 34, y + 11), FEATURE_LABELS.get(feature, feature), anchor="lm", fill=(0, 0, 0), font=font("regular", 17))

    out.parent.mkdir(parents=True, exist_ok=True)
    img.save(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    rows, features = collect(args.model_dir)
    if not rows:
        raise SystemExit(f"No ptFine cent7 BDT JSON/metadata pairs found in {args.model_dir}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    csv_path = args.outdir / "ptfine_cent7_bdt_feature_importance_by_route.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    draw_overall(rows, features, args.outdir / "ptfine_cent7_bdt_all_feature_gain_overall.png")
    draw_stacked_by_pt(rows, features, args.outdir / "ptfine_cent7_bdt_feature_gain_stacked_by_et.png")
    manifest = {
        "schema": "RJ_AUAU_BDT_FEATURE_IMPORTANCE_PLOTS_V1",
        "model_dir": str(args.model_dir),
        "n_routes": len({r["route"] for r in rows}),
        "n_features": len(features),
        "csv": str(csv_path),
        "plots": [
            str(args.outdir / "ptfine_cent7_bdt_all_feature_gain_overall.png"),
            str(args.outdir / "ptfine_cent7_bdt_feature_gain_stacked_by_et.png"),
        ],
    }
    (args.outdir / "ptfine_cent7_bdt_feature_importance_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
