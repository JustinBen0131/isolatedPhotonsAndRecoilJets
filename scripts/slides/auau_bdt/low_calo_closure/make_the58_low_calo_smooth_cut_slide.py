#!/usr/bin/env python3
"""THE-58 slide-5 (v3, valley-anchored): low-calo veto with ONE smooth T(c).

Fix vs v2: the cut is fit to the per-5%-bin VALLEY just above the fixed low-calo
artifact (log10 E_calo ~ 2.845), not to median-5*sigma. The artifact is pinned
at a fixed energy in every centrality (centrality-independent -> detector
artifact); the main band slides left toward it. median-5sigma got dragged left
by the sliding band and slid off the artifact above ~35% (under-cut 40-55%).

The valley anchor keeps a smooth continuous curve T(c) (a near-flat quadratic)
but tracks the artifact's right edge, so it removes ~6-7% in every bin where the
artifact is RESOLVED (0 to ~55%). Above ~55% the band merges with the artifact
and they cannot be separated by calo energy -> not cut there (shown combined).

Raw event counts (not normalized). x-axis every panel: log10 E_calo.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np

SRC = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
HIST_JSON = SRC / "the32_low_calo_full_count_histograms_v1.json"

OUTDIR = Path("dataOutput/auauTightBDTValidation/THE58_lowCaloSmoothCut_20260614")
PNG = OUTDIR / "the58_low_calo_smooth_cut_slide5_v3.png"
MANIFEST = OUTDIR / "the58_low_calo_smooth_cut_slide5_v3.manifest.json"

ARTIFACT = 2.845       # fixed low-calo artifact location in log10 E_calo
FIT_DEGREE = 2

BLUE = "#1F77B4"
RED = "#D62728"
INK = "#172033"
MUTED = "#5D6B7A"
BOX_EDGE = "#CBD5E1"
GOLD = "#B7791F"
GOLD_SOFT = "#FFF7ED"
GOLD_LINE = "#F1C27D"
PURPLE = "#7A5BBE"
ARTIFACT_COLOR = "#2E7E8E"  # dashed artifact line + the "fixed artifact" words (color-coded)

plt.rcParams.update(
    {"font.family": "Times New Roman", "mathtext.fontset": "stix", "axes.unicode_minus": False}
)


def step_xy(edges, counts):
    return np.repeat(edges, 2)[1:-1], np.repeat(counts, 2)


def comma(v):
    return f"{int(round(v)):,}"


def smooth(y, w=3):
    return np.convolve(y, np.ones(w) / w, mode="same")


def build():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    hist = json.loads(HIST_JSON.read_text())
    edges = np.asarray(hist["energy_edges"], dtype="float64")
    centers = 0.5 * (edges[:-1] + edges[1:])
    panels = hist["panels"]

    # --- per-bin valley anchor (right edge of the fixed artifact bump) ---
    cents, anchors, sepflags, totals = [], [], [], []
    for p in panels:
        total = np.asarray(p["retained_hist"], float) + np.asarray(p["removed_hist"], float)
        s = smooth(total, 3)
        c = 0.5 * (p["cent_lo"] + p["cent_hi"])
        win = [i for i in range(2, len(s) - 2)
               if 2.80 <= centers[i] <= 2.90 and s[i] >= s[i - 1] and s[i] >= s[i + 1] and s[i] > 20]
        fi = win[int(np.argmax([s[i] for i in win]))] if win else None
        mi = int(np.argmax(s))
        separable = fi is not None and centers[mi] - centers[fi] > 0.06
        if separable:
            j = fi
            while j < mi - 1 and not (s[j] <= s[j + 1] and s[j] < s[j - 1]):
                j += 1
            anchors.append(centers[j])
        else:
            anchors.append(None)
        cents.append(c); sepflags.append(separable); totals.append(total)
    cents = np.array(cents)

    sep = np.array(sepflags)
    coef = np.polyfit(cents[sep], np.array([a for a in anchors if a is not None]), FIT_DEGREE)

    # --- apply the smooth cut; merged bins are not cut (cannot separate) ---
    rebuilt = []
    for c, total, s_ok in zip(cents, totals, sepflags):
        thr = float(np.polyval(coef, c))
        if s_ok:
            below = centers < thr
            removed = float(total[below].sum())
            rh = np.where(below, total, 0.0); kh = np.where(below, 0.0, total)
        else:
            removed = 0.0; rh = np.zeros_like(total); kh = total
        rebuilt.append(dict(lo=int(round(c - 2.5)), hi=int(round(c + 2.5)), c=c, thr=thr,
                            sep=s_ok, total=float(total.sum()), removed=removed,
                            removed_hist=rh, retained_hist=kh, full_hist=total,
                            frac=removed / max(total.sum(), 1.0)))

    central = [r for r in rebuilt if r["sep"]]
    periph = [r for r in rebuilt if not r["sep"]]
    comb = dict(lo=periph[0]["lo"], hi=periph[-1]["hi"], nbins=len(periph),
                full_hist=sum(r["full_hist"] for r in periph),
                total=sum(r["total"] for r in periph))
    tot_all = sum(r["total"] for r in rebuilt)
    rem_all = sum(r["removed"] for r in rebuilt)
    cen_ymax = max(max(r["retained_hist"].max(), r["removed_hist"].max()) for r in central)
    comb_ymax = comb["full_hist"].max()

    # ---------------- figure ----------------
    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(3, 4, left=0.052, right=0.965, top=0.578, bottom=0.094,
                          hspace=1.22, wspace=0.135)

    def style(ax):
        ax.set_xlim(edges[0], edges[-1])
        ax.set_xticks([2.8, 3.0, 3.2, 3.4])
        ax.set_xlabel(r"$\log_{10}E_{\rm calo}$", fontsize=9.8, color=MUTED, labelpad=3.0)
        ax.tick_params(axis="both", labelsize=8.4, colors=INK, length=3.0, width=0.85)
        for sp in ax.spines.values():
            sp.set_linewidth(0.95); sp.set_color(INK)

    def draw(ax, h, color, alpha):
        # one CONNECTED filled step over the populated extent (first..last
        # nonzero bin); empty interior bins drop to the baseline but stay
        # connected -> coherent histogram, no dashed gaps. Clipping to the
        # nonzero extent keeps bump and band as separate shapes (no baseline
        # line strung across the empty gap between them).
        h = np.asarray(h, float)
        nz = np.nonzero(h > 0)[0]
        if nz.size == 0:
            return
        i0, i1 = int(nz[0]), int(nz[-1])
        e = edges[i0:i1 + 2]
        c = np.maximum(h[i0:i1 + 1], 0.9)
        x, y = step_xy(e, c)
        ax.plot(x, y, color=color, lw=1.9)
        ax.fill_between(x, y, 0.9, color=color, alpha=alpha)

    def cut_panel(ax, r):
        ax.set_facecolor("white"); ax.set_yscale("log")
        draw(ax, r["retained_hist"], BLUE, 0.06)
        draw(ax, r["removed_hist"], RED, 0.10)
        ax.axvline(r["thr"], color="white", lw=3.8, ls=(0, (3.0, 3.0)), zorder=4)
        ax.axvline(r["thr"], color=GOLD, lw=2.1, ls=(0, (3.0, 3.0)), zorder=5)
        ax.set_ylim(0.8, max(12.0, cen_ymax * 2.6))
        ax.set_title(f"{r['lo']}-{r['hi']}%", fontsize=11.4, fontweight="bold", color=INK, pad=15)
        ax.text(0.0, 1.05, rf"$N_{{\rm acc}}$={comma(r['total']-r['removed'])}", transform=ax.transAxes,
                fontsize=9.6, color=BLUE, fontweight="bold", ha="left", va="bottom", clip_on=False, zorder=10)
        ax.text(0.50, 1.05, rf"$N_{{\rm cut}}$={comma(r['removed'])}", transform=ax.transAxes,
                fontsize=9.6, color=RED, fontweight="bold", ha="left", va="bottom", clip_on=False, zorder=10)
        ax.text(1.0, 1.05, f"({100*r['frac']:.1f}%)", transform=ax.transAxes,
                fontsize=9.6, color=RED, fontweight="bold", ha="right", va="bottom", clip_on=False, zorder=10)
        style(ax)

    def merged_panel(ax, c):
        ax.set_facecolor("#FBF7FF"); ax.set_yscale("log")
        draw(ax, c["full_hist"], PURPLE, 0.07)
        ax.set_ylim(0.8, max(12.0, comb_ymax * 2.6))
        ax.set_title(f"{c['lo']}-{c['hi']}%  ({c['nbins']} bins)", fontsize=11.2,
                     fontweight="bold", color=INK, pad=15)
        ax.text(0.0, 1.05, rf"$N_{{\rm acc}}$={comma(c['total'])}", transform=ax.transAxes, fontsize=9.6,
                color=BLUE, fontweight="bold", ha="left", va="bottom", clip_on=False, zorder=10)
        box = FancyBboxPatch((0.505, 0.470), 0.465, 0.485,
                             boxstyle="round,pad=0.015,rounding_size=0.04",
                             transform=ax.transAxes, facecolor="#F3EDFB",
                             edgecolor=PURPLE, linewidth=2.0, zorder=8)
        ax.add_patch(box)
        ax.text(0.7375, 0.800, "no low energy", transform=ax.transAxes, fontsize=14.0,
                fontweight="bold", color=PURPLE, ha="center", va="center", zorder=9)
        ax.text(0.7375, 0.625, "events cut", transform=ax.transAxes, fontsize=14.0,
                fontweight="bold", color=PURPLE, ha="center", va="center", zorder=9)
        style(ax)

    slots = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3), (2, 0), (2, 1), (2, 2)]
    leftcol = []
    for (rr, cc), r in zip(slots, central):
        ax = fig.add_subplot(gs[rr, cc]); cut_panel(ax, r)
        if cc == 0:
            leftcol.append(ax)
    axm = fig.add_subplot(gs[2, 3]); merged_panel(axm, comb)
    for ax in leftcol:
        ax.set_ylabel("raw event counts", fontsize=8.8, color=INK)

    # ---------------- header ----------------
    fig.text(0.052, 0.966, "Low-calo event veto: smooth threshold vs centrality",
             fontsize=22.0, fontweight="bold", color=INK, va="top")
    fig.text(0.052, 0.922, "Raw counts, not normalized.", fontsize=11.8, color=MUTED, va="top")
    fig.text(0.226, 0.922,
             r"Dashed line = fixed low-calo artifact ($\log_{10}E_{\rm calo}\approx2.85$, every centrality).",
             fontsize=11.8, color=ARTIFACT_COLOR, fontweight="bold", va="top")

    chip = FancyBboxPatch((0.052, 0.842), 0.605, 0.044,
                          boxstyle="round,pad=0.004,rounding_size=0.006",
                          transform=fig.transFigure, facecolor=GOLD_SOFT,
                          edgecolor=GOLD_LINE, linewidth=1.0, zorder=1)
    fig.patches.append(chip)
    fig.text(0.062, 0.864, "Veto:", fontsize=13.0, fontweight="bold", color=GOLD, va="center", zorder=3)
    fig.text(0.111, 0.864,
             r"$\mathrm{remove\ if}\ \log_{10}E_{\rm calo}<T(c),\ \ E_{\rm calo}=E_{\rm CEMC}+E_{\rm IHCal}+E_{\rm OHCal}+1$",
             fontsize=12.2, color=INK, va="center", ha="left", zorder=3)

    # two large arrowhead bullets, bold lead phrase
    fig.text(0.058, 0.812, r"$\blacktriangleright$", fontsize=11.5, color=GOLD, va="center", zorder=3)
    fig.text(0.086, 0.812, "Removes ~6-7%", fontsize=14.0, fontweight="bold", color=INK, va="center", zorder=3)
    fig.text(0.205, 0.812, r"wherever the artifact is resolved ($T(c)$ tracks the per-bin valley).",
             fontsize=13.2, color=INK, va="center", zorder=3)
    fig.text(0.058, 0.762, r"$\blacktriangleright$", fontsize=11.5, color=GOLD, va="center", zorder=3)
    fig.text(0.086, 0.762, "Past ~55%:", fontsize=14.0, fontweight="bold", color=INK, va="center", zorder=3)
    fig.text(0.168, 0.762, "band merges with the artifact, so it cannot be separated - not cut.",
             fontsize=13.2, color=INK, va="center", zorder=3)

    # removed-before-training takeaway (soft, not blinding)
    rembox = FancyBboxPatch((0.052, 0.677), 0.404, 0.048,
                            boxstyle="round,pad=0.006,rounding_size=0.008",
                            transform=fig.transFigure, facecolor="white",
                            edgecolor="#C96A6A", linewidth=1.8, zorder=2)
    fig.patches.append(rembox)
    fig.text(0.066, 0.701, "Removed before training:", fontsize=14.0, fontweight="bold",
             color="#B23B3B", va="center", zorder=3)
    fig.text(0.286, 0.701, f"{comma(rem_all)} / {comma(tot_all)}  ({100*rem_all/tot_all:.1f}%)",
             fontsize=14.0, fontweight="bold", color=INK, va="center", zorder=3)

    # legend (large line samples, spread out)
    ly = 0.640
    for x0, col, dash, lab, tx in [
        (0.058, BLUE, None, "kept", 0.096),
        (0.190, RED, None, "removed", 0.228),
        (0.346, GOLD, (0, (4, 3)), r"$T(c)$ cut", 0.384),
    ]:
        ln = Line2D([x0, x0 + 0.032], [ly, ly], color=col, lw=3.6, transform=fig.transFigure)
        if dash:
            ln.set_linestyle(dash)
        fig.add_artist(ln)
        fig.text(tx, ly, lab, fontsize=13.0, color=INK, va="center")

    # ---- T(c) curve, top-right ----
    cbox = FancyBboxPatch((0.672, 0.690), 0.293, 0.270,
                          boxstyle="round,pad=0.006,rounding_size=0.008",
                          transform=fig.transFigure, facecolor="white",
                          edgecolor=BOX_EDGE, linewidth=1.1, zorder=1)
    fig.patches.append(cbox)
    fig.text(0.818, 0.948, "One smooth veto curve  T(c)", fontsize=11.4, fontweight="bold",
             color=INK, ha="center", va="top", zorder=3)
    ins = fig.add_axes([0.714, 0.804, 0.226, 0.104])
    ins.set_zorder(6); ins.set_facecolor("white")
    cc = np.linspace(0, 80, 200)
    tv = np.polyval(coef, cc)
    ins.axhline(ARTIFACT, color=ARTIFACT_COLOR, lw=1.6, ls=(0, (2, 2)), zorder=2)
    ins.text(2, ARTIFACT - 0.013, "fixed artifact", fontsize=7.8, fontweight="bold",
             color=ARTIFACT_COLOR, va="top", zorder=3)
    above = tv >= ARTIFACT
    ins.plot(cc[above], tv[above], color=GOLD, lw=2.7, zorder=4)
    ins.plot(cc[~above], tv[~above], color=GOLD, lw=2.0, alpha=0.30, zorder=4)
    ins.scatter(cents[sep], [a for a in anchors if a is not None], s=12, color=INK, zorder=5)
    ins.text(78, 2.972, "cut sinks below\nartifact: 0 removed", fontsize=6.8, color=MUTED,
             ha="right", va="top", linespacing=1.2, zorder=6)
    ins.set_xlim(0, 80); ins.set_ylim(2.73, 3.0)
    ins.set_xlabel("centrality (%)", fontsize=8.4, color=INK, labelpad=2.0)
    ins.set_ylabel(r"$T(c)$, log$_{10}E_{\rm calo}$", fontsize=8.2, color=INK, labelpad=2.0)
    ins.tick_params(labelsize=7.4, colors=INK, length=2.6, width=0.8)
    for sp in ins.spines.values():
        sp.set_linewidth(0.85); sp.set_color(BOX_EDGE)
    a2, a1, a0 = coef
    e2 = int(np.floor(np.log10(abs(a2)))); m2 = a2 / 10.0 ** e2
    e1 = int(np.floor(np.log10(abs(a1)))); m1 = a1 / 10.0 ** e1
    eqbox = FancyBboxPatch((0.682, 0.692), 0.272, 0.054,
                           boxstyle="round,pad=0.004,rounding_size=0.006",
                           transform=fig.transFigure, facecolor=GOLD_SOFT,
                           edgecolor=GOLD_LINE, linewidth=1.4, zorder=6)
    fig.patches.append(eqbox)
    eqn = (rf"$T(c)={m2:.1f}\times10^{{{e2}}}\,c^{{2}}"
           rf"{'+' if a1>=0 else '-'}{abs(m1):.1f}\times10^{{{e1}}}\,c"
           rf"{'+' if a0>=0 else '-'}{abs(a0):.2f}$")
    fig.text(0.818, 0.719, eqn, fontsize=15.0, color=INK, ha="center", va="center", zorder=8)

    fig.savefig(PNG, dpi=200, facecolor="white")
    plt.close(fig)
    MANIFEST.write_text(json.dumps(dict(
        png=str(PNG), method="valley_anchor_smooth_fit", artifact=ARTIFACT, degree=FIT_DEGREE,
        coef=[float(x) for x in coef], total=tot_all, removed=rem_all, removed_frac=rem_all / tot_all,
        per_bin=[dict(lo=r["lo"], hi=r["hi"], thr=float(r["thr"]), sep=bool(r["sep"]), frac=float(r["frac"])) for r in rebuilt],
    ), indent=2))
    print(f"wrote {PNG}")
    print(f"overall removed {comma(rem_all)}/{comma(tot_all)} = {100*rem_all/tot_all:.2f}%")
    print(f"separable (cut) bins: {len(central)}; merged bins: {len(periph)}")
    print(f"T(c) coef: {[float(x) for x in coef]}")


if __name__ == "__main__":
    build()
