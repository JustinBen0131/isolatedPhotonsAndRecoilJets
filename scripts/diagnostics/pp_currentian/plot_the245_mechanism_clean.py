#!/usr/bin/env python3
"""THE-245 single-panel result.

Cause on x, effect on y, for identical candidates under two topo configurations:

    x = centroid dR between the photon and the topo cluster containing its
        leading CEMC tower
    y = topo-cluster isolation actually assigned to that candidate

If the negative branch is produced by the containing cluster's centroid leaving
the isolation cone, isolation collapses exactly at x = 0.4 and nowhere else.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

INK = "#1f2933"
TEAL = "#2a7f8f"
RED = "#d1495b"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    df = df[df["resolved"] == 1]
    ppg12 = df[df["ppg12_n_containing"] > 0]
    sam = df[df["sam_n_containing"] > 0]

    fig, ax = plt.subplots(figsize=(9.6, 6.2))

    ax.axhline(0.0, lw=0.9, color=INK, alpha=0.35, zorder=1)
    ax.axvline(0.4, ls="--", lw=1.6, color=INK, alpha=0.8, zorder=1)

    ax.scatter(ppg12["ppg12_match_dr"], ppg12["ppg12_legacy_iso04"],
               s=26, c=TEAL, alpha=0.75, linewidths=0, zorder=3,
               label=f"PPG12 settings   (n={len(ppg12)})")
    ax.scatter(sam["sam_match_dr"], sam["sam_legacy_iso04"],
               s=26, c=RED, alpha=0.7, linewidths=0, zorder=2,
               label=f"TopoClusterReco settings   (n={len(sam)})")

    ax.set_xlabel(r"centroid $\Delta R$  between the photon and the topo cluster containing it",
                  fontsize=11.5)
    ax.set_ylabel(r"topo-cluster isolation  [GeV]", fontsize=11.5)
    ax.set_xlim(-0.06, 3.6)
    ax.set_ylim(-24, 22)
    ax.grid(alpha=0.15, lw=0.6)
    ax.legend(loc="lower left", fontsize=10.5, framealpha=0.94)

    ax.annotate(r"$R = 0.4$ isolation cone", xy=(0.4, 18.5), xytext=(0.62, 18.5),
                fontsize=10.5, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.1))

    ax.set_title("pp Photon+jet10, 1000 events, identical candidates in both configurations",
                 fontsize=12.5, color=INK, pad=12)

    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
