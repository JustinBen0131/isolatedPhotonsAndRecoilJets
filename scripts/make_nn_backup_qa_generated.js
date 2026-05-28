#!/usr/bin/env node
/*
 * Generate a clean three-panel neural-network QA backup graphic from compact
 * CSV/JSON metrics, avoiding screenshot-like embedded legacy plots.
 */

const fs = require("fs");
const path = require("path");
let sharp;
try {
  sharp = require("sharp");
} catch (_) {
  sharp = require("/Users/patsfan753/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp");
}

const repo = path.resolve(__dirname, "..");
const outdir = path.join(repo, "dataOutput/auauMLPDiagnosticPlots/backup_nn_qa_triptych_20260519");
const outpng = path.join(outdir, "nn_backup_qa_three_plot_story_generated.png");

const historyCsv = path.join(
  repo,
  "dataOutput/auauMLPDiagnosticPlots/aligned_stack_score_and_training_20260513/training_curves/mlp_training_history_combined.csv"
);
const heatmapCsv = path.join(
  repo,
  "dataOutput/auauMLPDiagnosticPlots/expert_validation_pack_full_uncapped_20260514_2038_final/pt_cent_heatmap_summary.csv"
);
const rankCsv = path.join(
  repo,
  "dataOutput/auauBDTMLPStackFullStat/stacked_bdt_mlp_targeted_fullstat_ptFine15to35_cent7_full_20260513_123934/stacked_sweep_rank_table.csv"
);

function csvRows(file) {
  const text = fs.readFileSync(file, "utf8").trim();
  const lines = text.split(/\r?\n/);
  const header = lines.shift().split(",");
  return lines.map((line) => {
    const vals = line.split(",");
    const row = {};
    header.forEach((key, i) => {
      const raw = vals[i] ?? "";
      const num = Number(raw);
      row[key] = raw !== "" && Number.isFinite(num) ? num : raw;
    });
    return row;
  });
}

function esc(s) {
  return String(s)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;");
}

function text(x, y, body, size, opts = {}) {
  const weight = opts.weight ?? 400;
  const fill = opts.fill ?? "#0d1424";
  const anchor = opts.anchor ? ` text-anchor="${opts.anchor}"` : "";
  const italic = opts.italic ? " font-style=\"italic\"" : "";
  return `<text x="${x}" y="${y}" font-family="Times New Roman, Times, serif" font-size="${size}" font-weight="${weight}" fill="${fill}"${anchor}${italic}>${esc(body)}</text>`;
}

function linePath(points) {
  return points
    .map((p, i) => `${i === 0 ? "M" : "L"} ${p[0].toFixed(1)} ${p[1].toFixed(1)}`)
    .join(" ");
}

function labelNum(v) {
  if (Math.abs(v - Math.round(v)) < 1e-9) return String(Math.round(v));
  return String(v).replace(/\.0$/, "");
}

function binLabel(lo, hi) {
  return `${labelNum(lo)}-${labelNum(hi)}`;
}

function interpColor(v, lo, hi, stops) {
  const t = Math.max(0, Math.min(1, (v - lo) / (hi - lo || 1)));
  const scaled = t * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(scaled));
  const f = scaled - i;
  const a = stops[i], b = stops[i + 1];
  const rgb = a.map((x, j) => Math.round(x + f * (b[j] - x)));
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
}

function luminance(rgbText) {
  const m = rgbText.match(/\d+/g).map(Number);
  return 0.2126 * m[0] + 0.7152 * m[1] + 0.0722 * m[2];
}

function drawLossPanel(history, x, y, w, h) {
  const rows = history.filter((r) => r.label === "MLP stack input: deep primary ratios");
  const epochs = rows.map((r) => r.epoch);
  const train = rows.map((r) => r.train_loss);
  const val = rows.map((r) => r.validation_loss);
  const best = rows.reduce((a, b) => (b.validation_loss < a.validation_loss ? b : a), rows[0]);
  const xmin = Math.min(...epochs), xmax = Math.max(...epochs);
  const ymin = Math.min(...train, ...val) - 0.003;
  const ymax = Math.max(...train, ...val) + 0.003;
  const px = x + 58, py = y + 135, pw = w - 98, ph = 430;
  const X = (e) => px + ((e - xmin) / (xmax - xmin)) * pw;
  const Y = (v) => py + ph - ((v - ymin) / (ymax - ymin)) * ph;
  const trainPts = rows.map((r) => [X(r.epoch), Y(r.train_loss)]);
  const valPts = rows.map((r) => [X(r.epoch), Y(r.validation_loss)]);
  let s = "";
  s += text(x + 24, y + 44, "1. Training stability", 30, { weight: 700 });
  s += text(x + 24, y + 78, "Validation curve plateaus while training loss keeps falling.", 19, { fill: "#565d6a" });
  s += `<rect x="${px}" y="${py}" width="${pw}" height="${ph}" fill="white" stroke="#c9d1dc"/>`;
  for (let i = 0; i <= 4; i += 1) {
    const yy = py + (i / 4) * ph;
    s += `<line x1="${px}" x2="${px + pw}" y1="${yy}" y2="${yy}" stroke="#e6e9ef"/>`;
    const valTick = ymax - (i / 4) * (ymax - ymin);
    s += text(px - 10, yy + 5, valTick.toFixed(3), 15, { fill: "#4c5565", anchor: "end" });
  }
  for (let e = 0; e <= 160; e += 40) {
    const xx = X(e || xmin);
    s += `<line x1="${xx}" x2="${xx}" y1="${py}" y2="${py + ph}" stroke="#eef1f5"/>`;
    s += text(xx, py + ph + 26, e, 16, { fill: "#4c5565", anchor: "middle" });
  }
  s += `<path d="${linePath(trainPts)}" fill="none" stroke="#1976b9" stroke-width="3"/>`;
  s += `<path d="${linePath(valPts)}" fill="none" stroke="#f28e2b" stroke-width="3"/>`;
  s += `<line x1="${X(best.epoch)}" x2="${X(best.epoch)}" y1="${py}" y2="${py + ph}" stroke="#6b7280" stroke-dasharray="6 5"/>`;
  s += text(X(best.epoch) + 8, py + 24, `best val epoch ${best.epoch}`, 16, { fill: "#4b5563" });
  s += `<line x1="${px + 18}" x2="${px + 54}" y1="${py + 28}" y2="${py + 28}" stroke="#1976b9" stroke-width="4"/>`;
  s += text(px + 64, py + 33, "training loss", 17, { fill: "#1f2937" });
  s += `<line x1="${px + 18}" x2="${px + 54}" y1="${py + 56}" y2="${py + 56}" stroke="#f28e2b" stroke-width="4"/>`;
  s += text(px + 64, py + 61, "validation loss", 17, { fill: "#1f2937" });
  s += text(px + pw / 2, py + ph + 58, "Epoch", 18, { fill: "#1f2937", anchor: "middle" });
  s += text(x + 24, y + h - 56, `Readout: validation loss reaches ${best.validation_loss.toFixed(3)} and then changes slowly.`, 20, { fill: "#0d1424", weight: 700 });
  return s;
}

function drawHeatmapPanel(rows, metric, x, y, w, h, titleStr, subtitle, lo, hi, stops, suffix = "") {
  const stack = rows.filter((r) => r.model === "BDT+MLP NN stack");
  const etBins = [...new Map(stack.map((r) => [`${r.et_lo}-${r.et_hi}`, [r.et_lo, r.et_hi]])).values()]
    .sort((a, b) => a[0] - b[0]);
  const centBins = [...new Map(stack.map((r) => [`${r.cent_lo}-${r.cent_hi}`, [r.cent_lo, r.cent_hi]])).values()]
    .sort((a, b) => a[0] - b[0]);
  const cellW = Math.floor((w - 120) / etBins.length);
  const cellH = Math.floor((h - 190) / centBins.length);
  const hx = x + 78;
  const hy = y + 122;
  let s = "";
  s += text(x + 24, y + 44, titleStr, 30, { weight: 700 });
  s += text(x + 24, y + 78, subtitle, 19, { fill: "#565d6a" });
  for (let iy = 0; iy < centBins.length; iy += 1) {
    const [clo, chi] = centBins[iy];
    s += text(hx - 12, hy + iy * cellH + cellH / 2 + 6, binLabel(clo, chi), 15, { fill: "#374151", anchor: "end" });
    for (let ix = 0; ix < etBins.length; ix += 1) {
      const [elo, ehi] = etBins[ix];
      const row = stack.find((r) => r.et_lo === elo && r.et_hi === ehi && r.cent_lo === clo && r.cent_hi === chi);
      const v = row ? row[metric] : NaN;
      const fill = Number.isFinite(v) ? interpColor(v, lo, hi, stops) : "#e5e7eb";
      const tx = hx + ix * cellW;
      const ty = hy + iy * cellH;
      s += `<rect x="${tx}" y="${ty}" width="${cellW - 2}" height="${cellH - 2}" fill="${fill}" stroke="white" stroke-width="1"/>`;
      const tfill = luminance(fill) < 120 ? "#ffffff" : "#111827";
      s += text(tx + cellW / 2, ty + cellH / 2 + 6, Number.isFinite(v) ? `${v.toFixed(2)}${suffix}` : "", 16, { fill: tfill, anchor: "middle", weight: 700 });
    }
  }
  etBins.forEach(([elo, ehi], ix) => {
    s += text(hx + ix * cellW + cellW / 2, hy + centBins.length * cellH + 26, binLabel(elo, ehi), 15, { fill: "#374151", anchor: "middle" });
  });
  s += text(hx + (etBins.length * cellW) / 2, hy + centBins.length * cellH + 55, "Candidate ET bin [GeV]", 17, { fill: "#1f2937", anchor: "middle" });
  return s;
}

async function main() {
  fs.mkdirSync(outdir, { recursive: true });
  const history = csvRows(historyCsv);
  const heatRows = csvRows(heatmapCsv);
  const rankRows = csvRows(rankCsv);
  const stackAll = rankRows.find((r) => r.model === "ptFine15to35_cent7_full_nn" && r.split === "all");

  const W = 1920, H = 1080;
  let svg = `<svg width="${W}" height="${H}" xmlns="http://www.w3.org/2000/svg">
  <rect width="${W}" height="${H}" fill="white"/>
  ${text(52, 68, "Neural-network QA backup: generated from validation metrics", 48, { weight: 700 })}
  ${text(54, 116, `BDT+MLP NN-stack diagnostic: AUC ${stackAll.auc.toFixed(3)}, WP80 fake ${stackAll.wp80_fake.toFixed(3)}; each panel answers a common review question.`, 25, { fill: "#565d6a" })}
  <rect x="50" y="154" width="555" height="850" rx="14" fill="#fafbfd" stroke="#d9dee7" stroke-width="2"/>
  <rect x="635" y="154" width="610" height="850" rx="14" fill="#fafbfd" stroke="#d9dee7" stroke-width="2"/>
  <rect x="1275" y="154" width="595" height="850" rx="14" fill="#fafbfd" stroke="#d9dee7" stroke-width="2"/>
  ${drawLossPanel(history, 50, 154, 555, 850)}
  ${drawHeatmapPanel(
    heatRows,
    "auc",
    635,
    154,
    610,
    850,
    "2. AUC stability",
    "Fine ET x centrality check for phase-space weak spots.",
    0.74,
    0.91,
    [
      [50, 64, 145],
      [53, 132, 187],
      [41, 176, 139],
      [170, 220, 70],
      [248, 221, 74],
    ]
  )}
  ${drawHeatmapPanel(
    heatRows,
    "wp80_fake",
    1275,
    154,
    595,
    850,
    "3. WP80 fake rate",
    "Fixed 80% signal efficiency: where backgrounds survive.",
    0.15,
    0.58,
    [
      [30, 64, 120],
      [82, 101, 138],
      [142, 133, 103],
      [214, 185, 65],
      [248, 221, 74],
    ]
  )}
  ${text(54, 1044, "Source: full-stat validation caches; stacker is backup/diagnostic until runtime parity and closure are finalized.", 18, { fill: "#666c77" })}
</svg>`;

  await sharp(Buffer.from(svg)).png().toFile(outpng);
  console.log(outpng);
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
