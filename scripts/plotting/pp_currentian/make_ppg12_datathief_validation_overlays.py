#!/usr/bin/env python3
"""Build PPG12 DataThief validation overlays.

This helper uses the official DataThief Java jar as the coordinate transform and
export engine. Pixel centers are detected from preserved PPG12 figure images,
then passed through DataThief with explicit axis reference points.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DATATHIEF_JAR = REPO / "local_tools/datathief/Datathief.jar"

FIG13_DIR = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/"
    / "shower_shape_reference_validation/fig13_e11_e33"
)
FIG13_IMAGE = FIG13_DIR / "ian_v4_fig13_e11_to_e33_crop.png"
FIG13_SDCC = FIG13_DIR / "ppg12_sdcc_fig13_e11_to_e33_histograms.json"
FIG13_PRIOR_DIGITIZER_SUMMARY = FIG13_DIR / "ppg12_ian_png_digitized_black_data_vs_sdcc_summary.json"
FIG13_OUT = FIG13_DIR / "ppg12_ian_datathief_export_black_data_vs_sdcc_overlay_slidefit_772x998.png"
FIG13_CSV = FIG13_DIR / "ppg12_ian_datathief_export_black_data_vs_sdcc_points.csv"
FIG13_RAW_CSV = FIG13_DIR / "ppg12_ian_datathief_export_black_data_vs_sdcc_points_raw_y18.csv"
FIG13_MANIFEST = FIG13_DIR / "ppg12_ian_datathief_export_black_data_vs_sdcc_manifest.json"

EFF_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/efficiency_stage"
FIG6_IMAGE = Path(
    "/var/folders/l3/f02nw86n5cn0tpf_zstf0ypr0000gn/T/TemporaryItems/"
    "NSIRD_screencaptureui_ZbDYpm/Screenshot 2026-06-30 at 11.17.20\u202fAM.png"
)
FIG6_SOURCE_COPY = EFF_DIR / "ppg12_fig6_paper_source_for_datathief_plotonly.png"
FIG6_SDCC = EFF_DIR / "ppg12_bdt_nom_efficiency_tefficiency_readback_full_remote_clean.csv"
FIG6_OUT = EFF_DIR / "ppg12_fig6_datathief_export_vs_sdcc_overlay_ratio.png"
FIG6_CSV = EFF_DIR / "ppg12_fig6_datathief_export_vs_sdcc_points.csv"
FIG6_MANIFEST = EFF_DIR / "ppg12_fig6_datathief_export_vs_sdcc_manifest.json"

WORK_DIR = EFF_DIR / "datathief_work"


JAVA_SOURCE = r"""
import java.io.*;
import java.util.*;
import java.lang.reflect.*;

public class DatathiefPointExporter {
    static class Ref {
        int index;
        double px;
        double py;
        String vx;
        String vy;
    }

    static class Pt {
        String series;
        double px;
        double py;
    }

    static String[] splitCsv(String line) {
        return line.split(",", -1);
    }

    static List<Ref> readRefs(String path) throws Exception {
        ArrayList<Ref> refs = new ArrayList<Ref>();
        BufferedReader br = new BufferedReader(new FileReader(path));
        String line = br.readLine();
        while ((line = br.readLine()) != null) {
            if (line.trim().length() == 0) continue;
            String[] f = splitCsv(line);
            Ref r = new Ref();
            r.index = Integer.parseInt(f[0]);
            r.px = Double.parseDouble(f[1]);
            r.py = Double.parseDouble(f[2]);
            r.vx = f[3];
            r.vy = f[4];
            refs.add(r);
        }
        br.close();
        return refs;
    }

    static LinkedHashMap<String, ArrayList<Pt>> readPoints(String path) throws Exception {
        LinkedHashMap<String, ArrayList<Pt>> bySeries = new LinkedHashMap<String, ArrayList<Pt>>();
        BufferedReader br = new BufferedReader(new FileReader(path));
        String line = br.readLine();
        while ((line = br.readLine()) != null) {
            if (line.trim().length() == 0) continue;
            String[] f = splitCsv(line);
            Pt p = new Pt();
            p.series = f[0];
            p.px = Double.parseDouble(f[1]);
            p.py = Double.parseDouble(f[2]);
            if (!bySeries.containsKey(p.series)) bySeries.put(p.series, new ArrayList<Pt>());
            bySeries.get(p.series).add(p);
        }
        br.close();
        return bySeries;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 4) {
            System.err.println("Usage: DatathiefPointExporter <image> <refs.csv> <points.csv> <out.csv>");
            System.exit(2);
        }
        String imagePath = args[0];
        String refsPath = args[1];
        String pointsPath = args[2];
        String outPath = args[3];

        DtFrame frame = new DtFrame();
        frame.open(imagePath);
        for (int i = 0; i < 200 && !frame.canvas.loaded(); i++) {
            Thread.sleep(25);
        }
        frame.setAxes("00 lin X - lin Y");
        Field refField = DataCanvas.class.getDeclaredField("refPoint");
        refField.setAccessible(true);
        Coordinate[] refPoints = (Coordinate[]) refField.get(frame.canvas);
        for (Ref r : readRefs(refsPath)) {
            RefData rd = frame.canvas.getRefData(r.index);
            Coordinate coord = new Coordinate(r.px, r.py);
            rd.setPoint(coord);
            rd.setXdata(r.vx);
            rd.setYdata(r.vy);
            refPoints[r.index] = coord;
        }
        refField.set(frame.canvas, refPoints);
        Field refsField = DataCanvas.class.getDeclaredField("refs");
        refsField.setAccessible(true);
        RefData[] refs = (RefData[]) refsField.get(frame.canvas);
        Computer computer = new Computer(
            refPoints,
            false,
            refs,
            frame.canvas.getAxes(),
            frame.canvas.translators[0],
            frame.canvas.translators[1]
        );

        LinkedHashMap<String, ArrayList<Pt>> bySeries = readPoints(pointsPath);
        PrintWriter out = new PrintWriter(new FileWriter(outPath));
        out.println("series,point_index,datathief_x,datathief_y");
        for (Map.Entry<String, ArrayList<Pt>> entry : bySeries.entrySet()) {
            ArrayList<Pt> pts = entry.getValue();
            for (int i = 0; i < pts.size(); i++) {
                Tuple in = new Tuple(pts.get(i).px, pts.get(i).py);
                Tuple transformed = new Tuple();
                computer.compute(in, transformed);
                out.println(
                    entry.getKey() + "," + i + "," +
                    Double.toString(transformed.x) + "," +
                    Double.toString(transformed.y)
                );
            }
        }
        out.close();
        frame.dispose();
        System.exit(0);
    }
}
"""


def jar_md5() -> str:
    return hashlib.md5(DATATHIEF_JAR.read_bytes()).hexdigest()


def ensure_datathief_ready() -> None:
    if not DATATHIEF_JAR.exists():
        raise FileNotFoundError(f"Missing DataThief jar: {DATATHIEF_JAR}")
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    java_file = WORK_DIR / "DatathiefPointExporter.java"
    java_file.write_text(JAVA_SOURCE)
    subprocess.run(
        ["javac", "-cp", str(DATATHIEF_JAR), str(java_file)],
        check=True,
        cwd=str(REPO),
    )


def run_datathief_export(
    image: Path,
    refs: list[tuple[int, float, float, float, float]],
    points: list[tuple[str, float, float]],
    output_csv: Path,
) -> None:
    ensure_datathief_ready()
    refs_csv = output_csv.with_suffix(".refs.csv")
    points_csv = output_csv.with_suffix(".pixels.csv")
    refs_csv.parent.mkdir(parents=True, exist_ok=True)
    with refs_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ref_index", "pixel_x", "pixel_y", "value_x", "value_y"])
        for row in refs:
            writer.writerow(row)
    with points_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["series", "pixel_x", "pixel_y"])
        for row in points:
            writer.writerow(row)
    classpath = f"{DATATHIEF_JAR}:{WORK_DIR}"
    subprocess.run(
        [
            "java",
            "-cp",
            classpath,
            "DatathiefPointExporter",
            str(image),
            str(refs_csv),
            str(points_csv),
            str(output_csv),
        ],
        check=True,
        cwd=str(REPO),
        timeout=20,
    )


def component_centers(mask: np.ndarray) -> list[dict[str, float]]:
    h, w = mask.shape
    seen = np.zeros_like(mask, dtype=bool)
    out: list[dict[str, float]] = []
    for y in range(h):
        xs = np.where(mask[y] & (~seen[y]))[0]
        for x in xs:
            if seen[y, x] or not mask[y, x]:
                continue
            stack = [(y, int(x))]
            seen[y, x] = True
            pts: list[tuple[int, int]] = []
            while stack:
                cy, cx = stack.pop()
                pts.append((cy, cx))
                for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ny = cy + dy
                    nx = cx + dx
                    if 0 <= ny < h and 0 <= nx < w and mask[ny, nx] and not seen[ny, nx]:
                        seen[ny, nx] = True
                        stack.append((ny, nx))
            if len(pts) < 4:
                continue
            yy = np.asarray([p[0] for p in pts], dtype=float)
            xx = np.asarray([p[1] for p in pts], dtype=float)
            out.append(
                {
                    "n": float(len(pts)),
                    "x_mean": float(xx.mean()),
                    "y_mean": float(yy.mean()),
                    "x_min": float(xx.min()),
                    "x_max": float(xx.max()),
                    "y_min": float(yy.min()),
                    "y_max": float(yy.max()),
                }
            )
    return out


def read_fig13_sdcc() -> dict[str, np.ndarray]:
    with FIG13_SDCC.open() as f:
        data = json.load(f)["data"]
    return {
        "centers": np.asarray(data["centers"], dtype=float),
        "values": np.asarray(data["values"], dtype=float),
        "errors": np.asarray(data["errors"], dtype=float),
    }


def find_fig13_points(image: Path, centers: np.ndarray) -> tuple[dict[str, float], list[tuple[str, float, float]]]:
    if FIG13_PRIOR_DIGITIZER_SUMMARY.exists():
        with FIG13_PRIOR_DIGITIZER_SUMMARY.open() as f:
            prior = json.load(f)
        axis = prior["axis_calibration"]
        rows = prior["rows"]
        return (
            {
                "x_left": float(axis["x_left"]),
                "x_right": float(axis["x_right"]),
                "y_bottom_0": float(axis["y_at_zero"]),
                "y_at_018": float(axis["y_at_018"]),
                "y_top_frame": float(axis["y_top_frame"]),
                "source_image": prior.get("source_image", str(image)),
                "pixel_source": str(FIG13_PRIOR_DIGITIZER_SUMMARY),
            },
            [("data", float(row["pixel_x"]), float(row["pixel_y"])) for row in rows],
        )

    arr = np.asarray(Image.open(image).convert("RGB"))
    gray = np.dot(arr[..., :3], [0.299, 0.587, 0.114])
    dark = gray < 65
    col_counts = dark.sum(axis=0)
    row_counts = dark.sum(axis=1)
    x_candidates = np.where(col_counts > 0.9 * dark.shape[0])[0]
    y_candidates = np.where(row_counts > 0.65 * dark.shape[1])[0]
    x_left = float(x_candidates[0])
    x_right = float(x_candidates[-1])
    y_top_frame = float(y_candidates[0])
    y_bottom = float(y_candidates[1])

    # The preserved Fig. 13 crop has a clear 0.18 major tick just below the top
    # frame. Use that tick rather than the frame top.
    tick_counts = dark[:, int(x_left) : int(x_left) + 30].sum(axis=1)
    tick_rows = np.where(
        (tick_counts > 18)
        & (np.arange(len(tick_counts)) > y_top_frame)
        & (np.arange(len(tick_counts)) < y_bottom)
    )[0]
    y_at_018 = float(tick_rows[0])

    points: list[tuple[str, float, float]] = []
    for c in centers:
        x_pred = x_left + c * (x_right - x_left)
        xlo = max(int(x_left + 4), int(round(x_pred)) - 13)
        xhi = min(int(x_right - 4), int(round(x_pred)) + 13)
        ylo = int(max(y_at_018 + 58, 70))
        yhi = int(y_bottom - 1)
        sub = dark[ylo : yhi + 1, xlo : xhi + 1]

        best: tuple[float, float, float] | None = None
        radius = 8
        yy, xx = np.ogrid[-radius : radius + 1, -radius : radius + 1]
        disk = xx * xx + yy * yy <= radius * radius
        for cy in range(ylo + radius, yhi - radius + 1):
            for cx in range(xlo + radius, xhi - radius + 1):
                patch = dark[cy - radius : cy + radius + 1, cx - radius : cx + radius + 1]
                density = float(patch[disk].mean())
                central = float(dark[cy - 3 : cy + 4, cx - 3 : cx + 4].mean())
                score = density + central - 0.005 * abs(cx - x_pred)
                if best is None or score > best[0]:
                    best = (score, float(cx), float(cy))
        if best and best[0] > 1.1:
            points.append(("data", best[1], best[2]))
            continue

        # Tail markers are clipped against the axis in the screenshot.
        yy_idx, xx_idx = np.where(sub)
        if len(xx_idx) == 0:
            raise RuntimeError(f"Could not find Fig. 13 marker near x={c}")
        points.append(("data", float(xlo + np.median(xx_idx)), float(ylo + np.median(yy_idx))))

    cal = {
        "x_left": x_left,
        "x_right": x_right,
        "y_bottom_0": y_bottom,
        "y_at_018": y_at_018,
        "y_top_frame": y_top_frame,
        "source_image": str(image),
    }
    return cal, points


def load_datathief_csv(path: Path) -> dict[str, list[tuple[float, float]]]:
    out: dict[str, list[tuple[float, float]]] = defaultdict(list)
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            out[row["series"]].append((float(row["datathief_x"]), float(row["datathief_y"])))
    return out


def write_fig13_outputs() -> None:
    sdcc = read_fig13_sdcc()
    cal, points = find_fig13_points(FIG13_IMAGE, sdcc["centers"])
    # DataThief's expression parser/export path rounds tiny axis ranges too
    # aggressively here, so use a 100x y-axis inside DataThief and divide the
    # exported y values back to the plotted 0-0.18 range below.
    datathief_y_scale = 100.0
    refs = [
        (0, cal["x_left"], cal["y_bottom_0"], 0.0, 0.0),
        (1, cal["x_right"], cal["y_bottom_0"], 1.0, 0.0),
        (2, cal["x_left"], cal["y_at_018"], 0.0, 0.18 * datathief_y_scale),
    ]
    run_datathief_export(FIG13_IMAGE, refs, points, FIG13_RAW_CSV)
    dt = load_datathief_csv(FIG13_RAW_CSV)["data"]
    dt_x = np.asarray([p[0] for p in dt], dtype=float)
    dt_y = np.asarray([p[1] / datathief_y_scale for p in dt], dtype=float)
    ratio = np.divide(dt_y, sdcc["values"], out=np.full_like(dt_y, np.nan), where=sdcc["values"] > 0)
    with FIG13_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["series", "point_index", "datathief_x", "datathief_y"])
        for i, (xv, yv) in enumerate(zip(dt_x, dt_y)):
            writer.writerow(["data", i, f"{xv:.12g}", f"{yv:.12g}"])

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.1,
        }
    )
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=200,
        sharex=True,
        gridspec_kw={"height_ratios": [3.25, 1.0], "hspace": 0.05},
    )
    ax.errorbar(
        sdcc["centers"],
        sdcc["values"],
        yerr=sdcc["errors"],
        fmt="o",
        color="black",
        ms=4.4,
        lw=1.0,
        label="SDCC ROOT data",
        zorder=3,
    )
    ax.plot(
        dt_x,
        dt_y,
        "s",
        color="#d62728",
        markerfacecolor="none",
        markeredgewidth=1.25,
        ms=4.8,
        label="DataThief export from IAN PNG",
        zorder=4,
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.18)
    ax.set_ylabel("normalized counts", fontsize=16)
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=12)
    ax.text(0.052, 0.92, "sPHENIX", transform=ax.transAxes, fontsize=14.5, fontstyle="italic", fontweight="bold")
    ax.text(0.235, 0.92, "Internal", transform=ax.transAxes, fontsize=14.5)
    ax.text(0.052, 0.845, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=11.5)
    ax.text(0.052, 0.785, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=11.5)
    ax.text(0.052, 0.725, r"$22<p_T<28$ GeV", transform=ax.transAxes, fontsize=11.5)
    ax.text(0.052, 0.665, "w/o nbkg cut", transform=ax.transAxes, fontsize=11.5)
    ax.legend(loc="upper right", frameon=False, fontsize=12.8, handlelength=1.3)

    rax.axhline(1.0, color="0.35", lw=1.0, ls=(0, (4, 4)))
    rax.plot(sdcc["centers"], ratio, "o", color="#d62728", ms=3.8)
    rax.set_ylim(0.88, 1.12)
    rax.set_ylabel("PNG / SDCC", fontsize=12.5)
    rax.set_xlabel(r"$e_{11}/e_{33}$", fontsize=15)
    rax.minorticks_on()
    rax.tick_params(which="both", direction="in", top=True, right=True, labelsize=11)
    fig.subplots_adjust(left=0.15, right=0.97, top=0.98, bottom=0.09)
    fig.savefig(FIG13_OUT)
    plt.close(fig)

    summary = {
        "artifact": str(FIG13_OUT),
        "source_image": str(FIG13_IMAGE),
        "pixel_source_summary": str(FIG13_PRIOR_DIGITIZER_SUMMARY) if FIG13_PRIOR_DIGITIZER_SUMMARY.exists() else None,
        "sdcc_json": str(FIG13_SDCC),
        "datathief_jar": str(DATATHIEF_JAR),
        "datathief_jar_md5": jar_md5(),
        "datathief_export_csv": str(FIG13_CSV),
        "datathief_raw_scaled_csv": str(FIG13_RAW_CSV),
        "datathief_internal_y_scale": datathief_y_scale,
        "axis_refs": refs,
        "axis_calibration_pixels": cal,
        "mean_abs_ratio_minus_one": float(np.nanmean(np.abs(ratio - 1.0))),
        "max_abs_ratio_minus_one": float(np.nanmax(np.abs(ratio - 1.0))),
        "note": "Marker pixel centers are from the saved Fig. 13 digitizer evidence; values are re-exported through the official DataThief jar coordinate transform.",
    }
    FIG13_MANIFEST.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def read_fig6_sdcc() -> dict[str, list[dict[str, float]]]:
    raw: dict[str, dict[float, dict[str, float]]] = defaultdict(dict)
    with FIG6_SDCC.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            stage = row["stage"]
            mid = float(row["pt_mid"])
            if mid < 10.0 or mid > 35.0:
                continue
            raw[stage][mid] = {
                "pt_low": float(row["pt_low"]),
                "pt_high": float(row["pt_high"]),
                "pt_mid": mid,
                "eff": float(row["efficiency"]),
                "err_low": float(row["err_low"]),
                "err_high": float(row["err_high"]),
            }
    out: dict[str, list[dict[str, float]]] = {"reco": [], "reco_id": [], "reco_id_iso": []}
    for mid in sorted(raw["reco"]):
        reco = raw["reco"][mid]
        ident = raw["id"][mid]
        all_eff = raw["all"][mid]
        reco_id_eff = reco["eff"] * ident["eff"]
        reco_id_err = reco_id_eff * math.hypot(
            ident["err_high"] / ident["eff"],
            reco["err_high"] / reco["eff"],
        )
        out["reco"].append(reco)
        out["reco_id"].append(
            {
                "pt_low": reco["pt_low"],
                "pt_high": reco["pt_high"],
                "pt_mid": mid,
                "eff": reco_id_eff,
                "err_low": reco_id_err,
                "err_high": reco_id_err,
            }
        )
        out["reco_id_iso"].append(all_eff)
    return out


def find_fig6_points(image: Path, sdcc: dict[str, list[dict[str, float]]]) -> tuple[dict[str, float], list[tuple[str, float, float]]]:
    arr = np.asarray(Image.open(image).convert("RGB"))
    # Copy the source image into dataOutput so the manifest is not tied to a
    # transient screenshot path.
    if image != FIG6_SOURCE_COPY:
        FIG6_SOURCE_COPY.write_bytes(Path(image).read_bytes())
    arr = np.asarray(Image.open(FIG6_SOURCE_COPY).convert("RGB"))

    gray = np.dot(arr[..., :3], [0.299, 0.587, 0.114])
    dark = gray < 60
    col_counts = dark.sum(axis=0)
    row_counts = dark.sum(axis=1)
    x_left = float(np.where(col_counts > 250)[0][0])
    x_frame_right = float(np.where(col_counts > 250)[0][-1])
    y_top = float(np.where(row_counts > 400)[0][0])
    y_bottom = float(np.where(row_counts > 400)[0][-1])

    # PPG12 Figure 6 has the 35 GeV major tick inside the frame. Use it for the
    # x-reference rather than the right frame edge.
    bottom_counts = dark[int(y_bottom) - 30 : int(y_bottom) + 1, :].sum(axis=0)
    long_ticks = np.where((bottom_counts > 15) & (np.arange(len(bottom_counts)) > x_left + 1) & (np.arange(len(bottom_counts)) < x_frame_right))[0]
    x_35 = float(long_ticks[-1])
    left_counts = dark[:, int(x_left) : int(x_left) + 35].sum(axis=1)
    major_y = np.where((left_counts > 17) & (np.arange(len(left_counts)) > y_top) & (np.arange(len(left_counts)) < y_bottom))[0]
    y_at_1 = float(major_y[0])

    x_ref_10 = x_left
    x_ref_35 = x_35
    stage_masks = {
        "reco": {
            "mask": dark,
            "search_y": (int(y_at_1) + 22, int(y_at_1) + 70),
            "x_half_width": 11,
        },
        "reco_id": {
            "mask": (arr[:, :, 0] > 150) & (arr[:, :, 1] < 130) & (arr[:, :, 2] > 95),
            "search_y": (245, 380),
            "x_half_width": 12,
        },
        "reco_id_iso": {
            "mask": (arr[:, :, 1] > 95) & (arr[:, :, 0] < 150) & (arr[:, :, 2] < 150),
            "search_y": (365, 510),
            "x_half_width": 12,
        },
    }
    points: list[tuple[str, float, float]] = []
    for stage, rows in sdcc.items():
        spec = stage_masks[stage]
        mask = spec["mask"]
        ylo, yhi = spec["search_y"]
        for row in rows:
            mid = row["pt_mid"]
            if mid < 10.0 or mid > 35.0:
                continue
            x_pred = x_ref_10 + (mid - 10.0) / 25.0 * (x_ref_35 - x_ref_10)
            xlo = max(0, int(round(x_pred)) - int(spec["x_half_width"]))
            xhi = min(mask.shape[1] - 1, int(round(x_pred)) + int(spec["x_half_width"]))
            sub = mask[ylo : yhi + 1, xlo : xhi + 1]
            yy, xx = np.where(sub)
            if len(xx) == 0:
                raise RuntimeError(f"No {stage} marker pixels near pT={mid}")
            # Median is robust against the horizontal error bar because the
            # marker and bar share y, while text fragments are excluded by y/x.
            px = float(xlo + np.median(xx))
            py = float(ylo + np.median(yy))
            points.append((stage, px, py))

    cal = {
        "x_ref_10": x_ref_10,
        "x_ref_35": x_ref_35,
        "x_frame_right": x_frame_right,
        "y_bottom_0": y_bottom,
        "y_at_1": y_at_1,
        "y_top_frame": y_top,
        "source_image": str(FIG6_SOURCE_COPY),
    }
    return cal, points


def write_fig6_outputs() -> None:
    sdcc = read_fig6_sdcc()
    cal, points = find_fig6_points(FIG6_IMAGE, sdcc)
    datathief_x_span = 25.0
    refs = [
        (0, cal["x_ref_10"], cal["y_bottom_0"], 0.0, 0.0),
        (1, cal["x_ref_35"], cal["y_bottom_0"], 1.0, 0.0),
        (2, cal["x_ref_10"], cal["y_at_1"], 0.0, 1.0),
    ]
    run_datathief_export(FIG6_SOURCE_COPY, refs, points, FIG6_CSV)
    dt_raw = load_datathief_csv(FIG6_CSV)
    dt = {
        stage: [(10.0 + x * datathief_x_span, y) for x, y in rows]
        for stage, rows in dt_raw.items()
    }
    with FIG6_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["series", "point_index", "datathief_x", "datathief_y"])
        for stage, rows in dt.items():
            for i, (xv, yv) in enumerate(rows):
                writer.writerow([stage, i, f"{xv:.12g}", f"{yv:.12g}"])

    colors = {
        "reco": "#111111",
        "reco_id": "#d62d91",
        "reco_id_iso": "#2f9638",
    }
    labels = {
        "reco": r"$\varepsilon_{\mathrm{reco}}$",
        "reco_id": r"$\varepsilon_{\mathrm{reco}}\times\varepsilon_{\mathrm{ID}}$",
        "reco_id_iso": r"$\varepsilon_{\mathrm{reco}}\times\varepsilon_{\mathrm{ID}}\times\varepsilon_{\mathrm{iso}}$",
    }

    ratio_summary = {}
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.15,
        }
    )
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=200,
        sharex=True,
        gridspec_kw={"height_ratios": [3.25, 1.0], "hspace": 0.05},
    )
    for stage in ("reco", "reco_id", "reco_id_iso"):
        rows = sdcc[stage]
        x = np.asarray([r["pt_mid"] for r in rows], dtype=float)
        y = np.asarray([r["eff"] for r in rows], dtype=float)
        xerr = np.asarray([[r["pt_mid"] - r["pt_low"] for r in rows], [r["pt_high"] - r["pt_mid"] for r in rows]], dtype=float)
        yerr = np.asarray([[r["err_low"] for r in rows], [r["err_high"] for r in rows]], dtype=float)
        d = np.asarray(dt[stage], dtype=float)
        dt_x = d[:, 0]
        dt_y = d[:, 1]
        ax.errorbar(
            x,
            y,
            xerr=xerr,
            yerr=yerr,
            fmt="o",
            ms=6.0,
            mfc="white",
            mec=colors[stage],
            mew=1.55,
            ecolor=colors[stage],
            elinewidth=1.0,
            capsize=0,
            alpha=0.95,
            zorder=5,
        )
        ax.plot(
            dt_x,
            dt_y,
            "o",
            color=colors[stage],
            markerfacecolor=colors[stage],
            markeredgewidth=0.9,
            ms=4.0,
            linestyle="None",
            zorder=6,
        )
        ratio = dt_y / y
        ratio_summary[stage] = {
            "mean_abs_ratio_minus_one": float(np.nanmean(np.abs(ratio - 1.0))),
            "max_abs_ratio_minus_one": float(np.nanmax(np.abs(ratio - 1.0))),
        }
        rax.plot(x, ratio, "o", color=colors[stage], markerfacecolor=colors[stage], markeredgewidth=0.8, ms=4.2, linestyle="None")

    ax.set_xlim(10.0, 35.0)
    ax.set_ylim(0.0, 1.15)
    ax.set_ylabel("Efficiency", fontsize=17)
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=12)
    ax.text(0.055, 0.93, "sPHENIX", transform=ax.transAxes, fontsize=15.0, fontstyle="italic", fontweight="bold")
    ax.text(0.255, 0.93, "Internal", transform=ax.transAxes, fontsize=15.0)
    ax.text(0.055, 0.865, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=12.5)
    ax.text(0.735, 0.93, "PYTHIA8", transform=ax.transAxes, fontsize=15.0)
    ax.text(0.735, 0.865, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=12.5)
    source_handles = [
        Line2D([0], [0], marker="o", color="0.15", markerfacecolor="0.15", markeredgecolor="0.15", lw=0, ms=6.0, label="DataThief export (filled)"),
        Line2D([0], [0], marker="o", color="0.15", markerfacecolor="white", markeredgecolor="0.15", markeredgewidth=1.45, lw=0, ms=6.0, label="SDCC ROOT (open)"),
    ]
    color_handles = [
        Line2D([0], [0], marker="o", color=colors["reco"], markerfacecolor=colors["reco"], lw=0, ms=5.8, label=labels["reco"]),
        Line2D([0], [0], marker="o", color=colors["reco_id"], markerfacecolor=colors["reco_id"], lw=0, ms=5.8, label=labels["reco_id"]),
        Line2D([0], [0], marker="o", color=colors["reco_id_iso"], markerfacecolor=colors["reco_id_iso"], lw=0, ms=5.8, label=labels["reco_id_iso"]),
    ]
    source_legend = ax.legend(
        handles=source_handles,
        loc="lower left",
        bbox_to_anchor=(0.02, 0.185),
        frameon=False,
        fontsize=12.1,
        handlelength=1.0,
        borderpad=0.2,
        labelspacing=0.55,
    )
    ax.add_artist(source_legend)
    ax.legend(
        handles=color_handles,
        loc="lower left",
        bbox_to_anchor=(0.02, 0.015),
        frameon=False,
        fontsize=12.1,
        handlelength=1.0,
        borderpad=0.2,
        labelspacing=0.55,
    )

    rax.axhline(1.0, color="0.35", lw=1.0, ls=(0, (4, 4)))
    rax.set_ylim(0.92, 1.08)
    rax.set_ylabel("PNG / SDCC", fontsize=12.5)
    rax.set_xlabel(r"$E_{\mathrm{T}}^{\gamma,\mathrm{truth}}$ [GeV]", fontsize=15)
    rax.minorticks_on()
    rax.tick_params(which="both", direction="in", top=True, right=True, labelsize=11)
    fig.subplots_adjust(left=0.14, right=0.97, top=0.98, bottom=0.09)
    fig.savefig(FIG6_OUT)
    plt.close(fig)

    manifest = {
        "artifact": str(FIG6_OUT),
        "source_image": str(FIG6_SOURCE_COPY),
        "sdcc_csv": str(FIG6_SDCC),
        "datathief_jar": str(DATATHIEF_JAR),
        "datathief_jar_md5": jar_md5(),
        "datathief_export_csv": str(FIG6_CSV),
        "datathief_internal_x_range": [0.0, 1.0],
        "datathief_exported_x_range": [10.0, 35.0],
        "axis_refs": refs,
        "axis_calibration_pixels": cal,
        "ratio_summary": ratio_summary,
        "note": "Marker centers detected from the preserved PPG12 Figure 6 PNG; values exported through the official DataThief jar.",
    }
    FIG6_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    write_fig13_outputs()
    write_fig6_outputs()
    print(FIG13_OUT)
    print(FIG6_OUT)


if __name__ == "__main__":
    main()
