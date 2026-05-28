#!/usr/bin/env node
/*
 * Expert NN QA backup graphic from compact validation tables.
 * This intentionally uses generated charts rather than screenshots of older
 * plots, so the information density is high while labels remain readable.
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
const outpng = path.join(outdir, "nn_backup_qa_deep_validations_2x3.png");

const files = {
  calibration: path.join(repo, "dataOutput/auauMLPDiagnosticPlots/expert_validation_pack_full_uncapped_20260514_2038_final/calibration_summary.csv"),
  rocScore: path.join(repo, "dataOutput/auauMLPDiagnosticPlots/expert_validation_pack_full_uncapped_20260514_2038_final/roc_score_summary.csv"),
  ptCent: path.join(repo, "dataOutput/auauMLPDiagnosticPlots/expert_validation_pack_full_uncapped_20260514_2038_final/pt_cent_heatmap_summary.csv"),
  correlation: path.join(repo, "dataOutput/auauMLPDiagnosticPlots/bdt_mlp_score_correlation_20260515/bdt_mlp_score_correlation_summary.csv"),
  wpCells: path.join(repo, "dataOutput/auauBDTMLPStackPromotion/bdt_mlp_stack_nn_wp80_uncapped_diagnostic_20260514_1915_nnstack_wp80_uncapped/stack_working_points_target80_cells.csv"),
};

function parseCsv(file) {
  const text = fs.readFileSync(file, "utf8").trim();
  const lines = text.split(/\r?\n/);
  const header = lines.shift().split(",");
  return lines.filter(Boolean).map((line) => {
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
  return String(s).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;");
}

function text(x, y, body, size, opts = {}) {
  const weight = opts.weight ?? 400;
  const fill = opts.fill ?? "#0d1424";
  const anchor = opts.anchor ? ` text-anchor="${opts.anchor}"` : "";
  return `<text x="${x}" y="${y}" font-family="Times New Roman, Times, serif" font-size="${size}" font-weight="${weight}" fill="${fill}"${anchor}>${esc(body)}</text>`;
}

function pathLine(points) {
  return points.map((p, i) => `${i ? "L" : "M"} ${p[0].toFixed(1)} ${p[1].toFixed(1)}`).join(" ");
}

function labelNum(v) {
  return Math.abs(v - Math.round(v)) < 1e-9 ? `${Math.round(v)}` : `${v}`.replace(/\.0$/, "");
}

function binLabel(lo, hi) {
  return `${labelNum(lo)}-${labelNum(hi)}`;
}

function color(v, lo, hi, stops) {
  const t = Math.max(0, Math.min(1, (v - lo) / (hi - lo || 1)));
  const x = t * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(x));
  const f = x - i;
  const rgb = stops[i].map((a, j) => Math.round(a + f * (stops[i + 1][j] - a)));
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
}

function lum(rgb) {
  const m = rgb.match(/\d+/g).map(Number);
  return 0.2126 * m[0] + 0.7152 * m[1] + 0.0722 * m[2];
}

function panel(x, y, w, h, title, subtitle) {
  return `<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="12" fill="#fafbfd" stroke="#d8dee8" stroke-width="2"/>
${text(x + 20, y + 36, title, 25, { weight: 700 })}
${text(x + 20, y + 62, subtitle, 16, { fill: "#586173" })}`;
}

function axes(x, y, w, h, xmin, xmax, ymin, ymax, xticks, yticks) {
  let s = `<rect x="${x}" y="${y}" width="${w}" height="${h}" fill="white" stroke="#c9d1dc"/>`;
  xticks.forEach((v) => {
    const xx = x + ((v - xmin) / (xmax - xmin)) * w;
    s += `<line x1="${xx}" x2="${xx}" y1="${y}" y2="${y + h}" stroke="#eef1f5"/>`;
    s += text(xx, y + h + 20, v.toFixed(v < 1 ? 1 : 0), 13, { fill: "#4b5563", anchor: "middle" });
  });
  yticks.forEach((v) => {
    const yy = y + h - ((v - ymin) / (ymax - ymin)) * h;
    s += `<line x1="${x}" x2="${x + w}" y1="${yy}" y2="${yy}" stroke="#eef1f5"/>`;
    s += text(x - 8, yy + 5, v.toFixed(2), 13, { fill: "#4b5563", anchor: "end" });
  });
  return s;
}

function drawCalibration(rows, x, y, w, h) {
  const px = x + 74, py = y + 94, pw = w - 106, ph = h - 150;
  const models = ["MLP", "BDT+MLP NN stack"];
  const cols = { MLP: "#2b83ba", "BDT+MLP NN stack": "#d95f02" };
  let s = panel(x, y, w, h, "Reliability calibration", "Predicted score vs observed truth fraction.");
  s += axes(px, py, pw, ph, 0, 1, 0, 1, [0, 0.5, 1], [0, 0.5, 1]);
  s += `<path d="M ${px} ${py + ph} L ${px + pw} ${py}" stroke="#9aa3af" stroke-dasharray="6 5" fill="none"/>`;
  models.forEach((m) => {
    const pts = rows.filter((r) => r.model === m).map((r) => [px + r.confidence * pw, py + ph - r.signal_fraction * ph]);
    s += `<path d="${pathLine(pts)}" stroke="${cols[m]}" stroke-width="3" fill="none"/>`;
    pts.forEach((p) => { s += `<circle cx="${p[0]}" cy="${p[1]}" r="4" fill="${cols[m]}"/>`; });
  });
  s += text(px + 12, py + 22, "perfect calibration", 13, { fill: "#6b7280" });
  s += text(px + pw - 118, py + ph - 48, "MLP ECE 0.365", 15, { fill: cols.MLP, weight: 700 });
  s += text(px + pw - 160, py + ph - 24, "Stack ECE 0.307", 15, { fill: cols["BDT+MLP NN stack"], weight: 700 });
  s += text(px + pw / 2, y + h - 18, "Predicted score", 14, { anchor: "middle", fill: "#1f2937" });
  s += `<text x="${x + 20}" y="${py + ph / 2}" transform="rotate(-90 ${x + 20} ${py + ph / 2})" font-family="Times New Roman, Times, serif" font-size="14" fill="#1f2937" text-anchor="middle">Observed truth fraction</text>`;
  return s;
}

function drawCentralityBars(rows, x, y, w, h) {
  const px = x + 58, py = y + 98, pw = w - 86, ph = h - 156;
  const models = ["BDT", "MLP", "BDT+MLP NN stack"];
  const cols = { BDT: "#4e79a7", MLP: "#59a14f", "BDT+MLP NN stack": "#e15759" };
  const cent = [...new Map(rows.map((r) => [`${r.cent_lo}-${r.cent_hi}`, [r.cent_lo, r.cent_hi]])).values()];
  let s = panel(x, y, w, h, "Centrality ROC check", "AUC by coarse centrality for three score definitions.");
  s += axes(px, py, pw, ph, 0, cent.length, 0.65, 0.84, [], [0.70, 0.76, 0.82]);
  const groupW = pw / cent.length;
  const barW = groupW / 4.2;
  cent.forEach(([clo, chi], i) => {
    const cx = px + i * groupW;
    models.forEach((m, j) => {
      const r = rows.find((q) => q.model === m && q.cent_lo === clo && q.cent_hi === chi);
      const bh = ((r.auc - 0.65) / (0.84 - 0.65)) * ph;
      const bx = cx + groupW * 0.18 + j * barW;
      s += `<rect x="${bx}" y="${py + ph - bh}" width="${barW - 2}" height="${bh}" fill="${cols[m]}"/>`;
      s += text(bx + barW / 2, py + ph - bh - 5, r.auc.toFixed(2), 11, { anchor: "middle", fill: "#111827" });
    });
    s += text(cx + groupW / 2, py + ph + 21, binLabel(clo, chi), 13, { anchor: "middle", fill: "#374151" });
  });
  models.forEach((m, j) => {
    const lx = px + 20 + j * 112;
    s += `<rect x="${lx}" y="${py - 26}" width="15" height="15" fill="${cols[m]}"/>`;
    s += text(lx + 21, py - 13, m.replace("BDT+MLP NN stack", "Stack"), 13, { fill: "#1f2937" });
  });
  return s;
}

function drawCorrelation(rows, x, y, w, h) {
  const all = rows.filter((r) => r.population === "all" && String(r.selection).includes("<Et<"));
  const agg = new Map();
  all.forEach((r) => {
    const m = String(r.selection).match(/(\d+)<Et<(\d+)/);
    if (!m) return;
    const key = `${m[1]}-${m[2]}`;
    if (!agg.has(key)) agg.set(key, []);
    agg.get(key).push(r.pearson_r);
  });
  const bins = [...agg.keys()].sort((a, b) => Number(a.split("-")[0]) - Number(b.split("-")[0]));
  const vals = bins.map((b) => agg.get(b).reduce((a, c) => a + c, 0) / agg.get(b).length);
  const px = x + 56, py = y + 96, pw = w - 86, ph = h - 150;
  let s = panel(x, y, w, h, "BDT/NN score agreement", "Mean Pearson correlation by ET bin; disagreement exposes new information.");
  s += axes(px, py, pw, ph, 0, bins.length - 1, 0.45, 0.95, [], [0.5, 0.7, 0.9]);
  const pts = vals.map((v, i) => [px + (i / (bins.length - 1)) * pw, py + ph - ((v - 0.45) / 0.5) * ph]);
  s += `<path d="${pathLine(pts)}" stroke="#7b3294" stroke-width="4" fill="none"/>`;
  pts.forEach((p, i) => {
    s += `<circle cx="${p[0]}" cy="${p[1]}" r="5" fill="#7b3294"/>`;
    s += text(p[0], py + ph + 20, bins[i], 12, { anchor: "middle", fill: "#374151" });
  });
  s += text(px + pw / 2, y + h - 18, "Candidate ET bin [GeV]", 14, { anchor: "middle", fill: "#1f2937" });
  s += text(px + 10, py + 22, "lower = models disagree more", 14, { fill: "#4b5563" });
  return s;
}

function drawHeatmap(rows, x, y, w, h, cfg) {
  const etBins = [...new Map(rows.map((r) => [`${r[cfg.etLo]}-${r[cfg.etHi]}`, [r[cfg.etLo], r[cfg.etHi]]])).values()].sort((a, b) => a[0] - b[0]);
  const centBins = [...new Map(rows.map((r) => [`${r[cfg.centLo]}-${r[cfg.centHi]}`, [r[cfg.centLo], r[cfg.centHi]]])).values()].sort((a, b) => a[0] - b[0]);
  const px = x + 58, py = y + 91;
  const cellW = Math.floor((w - 92) / etBins.length);
  const cellH = Math.floor((h - 154) / centBins.length);
  let s = panel(x, y, w, h, cfg.title, cfg.subtitle);
  centBins.forEach(([clo, chi], iy) => {
    s += text(px - 9, py + iy * cellH + cellH / 2 + 5, binLabel(clo, chi), 12, { anchor: "end", fill: "#374151" });
    etBins.forEach(([elo, ehi], ix) => {
      const r = rows.find((q) => q[cfg.etLo] === elo && q[cfg.etHi] === ehi && q[cfg.centLo] === clo && q[cfg.centHi] === chi && (!cfg.model || q.model === cfg.model));
      const v = r?.[cfg.metric];
      const fill = Number.isFinite(v) ? color(v, cfg.lo, cfg.hi, cfg.stops) : "#e5e7eb";
      const tx = px + ix * cellW, ty = py + iy * cellH;
      s += `<rect x="${tx}" y="${ty}" width="${cellW - 2}" height="${cellH - 2}" fill="${fill}" stroke="white"/>`;
      s += text(tx + cellW / 2, ty + cellH / 2 + 5, Number.isFinite(v) ? v.toFixed(cfg.digits ?? 2) : "", 12, { anchor: "middle", fill: lum(fill) < 125 ? "#fff" : "#111827", weight: 700 });
    });
  });
  etBins.forEach(([elo, ehi], ix) => s += text(px + ix * cellW + cellW / 2, py + centBins.length * cellH + 19, binLabel(elo, ehi), 12, { anchor: "middle", fill: "#374151", weight: 700 }));
  s += text(px + (etBins.length * cellW) / 2, y + h - 12, "ET bin [GeV]", 13, { anchor: "middle", fill: "#1f2937" });
  return s;
}

async function main() {
  fs.mkdirSync(outdir, { recursive: true });
  const cal = parseCsv(files.calibration);
  const roc = parseCsv(files.rocScore);
  const ptCent = parseCsv(files.ptCent).filter((r) => r.model === "BDT+MLP NN stack");
  const corr = parseCsv(files.correlation);
  const wp = parseCsv(files.wpCells);

  const W = 1920, H = 1080;
  const panelW = 580, panelH = 390;
  const xs = [42, 670, 1298];
  const ys = [158, 590];

  let svg = `<svg width="${W}" height="${H}" xmlns="http://www.w3.org/2000/svg">
  <rect width="${W}" height="${H}" fill="white"/>
  ${text(48, 62, "Deep NN QA backup: calibration, disagreement, phase space, and WP stability", 43, { weight: 700 })}
  ${text(50, 108, "Generated from compact full-stat validation summaries; intended for detailed backup discussion, not the main story.", 24, { fill: "#586173" })}
  ${drawCalibration(cal, xs[0], ys[0], panelW, panelH)}
  ${drawCentralityBars(roc, xs[1], ys[0], panelW, panelH)}
  ${drawCorrelation(corr, xs[2], ys[0], panelW, panelH)}
  ${drawHeatmap(ptCent, xs[0], ys[1], panelW, panelH, {
    title: "Stack AUC map",
    subtitle: "Cell-by-cell ranking power after BDT+MLP stacking.",
    metric: "auc", model: "BDT+MLP NN stack",
    etLo: "et_lo", etHi: "et_hi", centLo: "cent_lo", centHi: "cent_hi",
    lo: 0.74, hi: 0.91, digits: 2,
    stops: [[50,64,145],[53,132,187],[41,176,139],[170,220,70],[248,221,74]],
  })}
  ${drawHeatmap(wp, xs[1], ys[1], panelW, panelH, {
    title: "WP80 threshold map",
    subtitle: "Runtime cut needed to keep 80% signal in each cell.",
    metric: "threshold",
    etLo: "pt_min", etHi: "pt_max", centLo: "centrality_min", centHi: "centrality_max",
    lo: 0.42, hi: 0.62, digits: 2,
    stops: [[57,81,162],[74,147,178],[69,181,137],[181,222,70],[253,231,88]],
  })}
  ${drawHeatmap(wp, xs[2], ys[1], panelW, panelH, {
    title: "WP80 fake-rate map",
    subtitle: "Background survival at the matched-efficiency working point.",
    metric: "background_fake_rate",
    etLo: "pt_min", etHi: "pt_max", centLo: "centrality_min", centHi: "centrality_max",
    lo: 0.12, hi: 0.47, digits: 2,
    stops: [[29,62,118],[82,101,138],[136,132,105],[201,174,73],[248,221,74]],
  })}
  ${text(50, 1042, "Caveat: stacker diagnostics use BDT score as an input and remain backup material until runtime parity and closure checks are complete.", 18, { fill: "#666c77" })}
  </svg>`;

  await sharp(Buffer.from(svg)).png().toFile(outpng);
  console.log(outpng);
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
