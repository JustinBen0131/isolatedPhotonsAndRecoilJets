#!/usr/bin/env node
/*
 * Compose the three highest-value NN QA plots into one backup display.
 * Uses the bundled Node/sharp runtime so the local analysis Python environment
 * does not need extra image packages.
 */

const path = require("path");
const fs = require("fs");
let sharp;
try {
  sharp = require("sharp");
} catch (_) {
  sharp = require("/Users/patsfan753/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp");
}

const repo = path.resolve(__dirname, "..");
const outdir = path.join(
  repo,
  "dataOutput/auauMLPDiagnosticPlots/backup_nn_qa_triptych_20260519"
);
const outpng = path.join(outdir, "nn_backup_qa_three_plot_story.png");
const outpngV2 = path.join(outdir, "nn_backup_qa_three_plot_story_wide_auc.png");

const inputs = [
  {
    path: path.join(
      repo,
      "dataOutput/auauMLPDiagnosticPlots/expert_validation_pack_full_uncapped_20260514_2038_final/01_training_validation_loss_mlp_and_nnstack.png"
    ),
    title: "1. Training stability",
    caption:
      "Train/validation loss checks whether later epochs improve generalization or only fit the training sample.",
    crop: { left: 70, top: 70, width: 2240, height: 1680 },
  },
  {
    path: path.join(
      repo,
      "dataOutput/auauMLPDiagnosticPlots/expert_validation_pack_full_uncapped_20260514_2038_final/07_auc_heatmap_pt_cent_mlp_vs_nnstack.png"
    ),
    title: "2. Phase-space performance",
    caption:
      "AUC by E_T and centrality shows where NN ranking is robust and where high-pT separation becomes harder.",
    crop: { left: 140, top: 90, width: 2860, height: 1080 },
  },
  {
    path: path.join(
      repo,
      "dataOutput/auauMLPDiagnosticPlots/expert_validation_pack_full_uncapped_20260514_2038_final/08_stack_wp80_fake_rate_grid.png"
    ),
    title: "3. WP80 operating behavior",
    caption:
      "At fixed 80% signal efficiency, the remaining background fake rate identifies bins that still drive purity risk.",
    crop: { left: 110, top: 70, width: 1290, height: 800 },
  },
];

function esc(s) {
  return s
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;");
}

function wrapWords(text, maxChars) {
  const words = text.split(/\s+/);
  const lines = [];
  let line = "";
  for (const word of words) {
    const candidate = line ? `${line} ${word}` : word;
    if (candidate.length <= maxChars || !line) {
      line = candidate;
    } else {
      lines.push(line);
      line = word;
    }
  }
  if (line) lines.push(line);
  return lines;
}

function textLines(lines, x, y, size, weight, color, lineHeight = 1.18) {
  return lines
    .map((line, idx) => {
      const dy = idx === 0 ? 0 : size * lineHeight;
      return `<text x="${x}" y="${y + dy}" font-family="Times New Roman, Times, serif" font-size="${size}" font-weight="${weight}" fill="${color}">${esc(line)}</text>`;
    })
    .join("\n");
}

async function trimBuffer(itemOrPath) {
  const imagePath = typeof itemOrPath === "string" ? itemOrPath : itemOrPath.path;
  const crop = typeof itemOrPath === "string" ? null : itemOrPath.crop;
  let pipeline = sharp(imagePath);
  if (crop) {
    pipeline = pipeline.extract(crop);
  }
  return pipeline
    .trim({ background: "#ffffff", threshold: 12 })
    .png()
    .toBuffer();
}

async function fitForPanel(imageBuffer, boxW, boxH) {
  const meta = await sharp(imageBuffer).metadata();
  const scale = Math.min(boxW / meta.width, boxH / meta.height);
  return {
    width: Math.round(meta.width * scale),
    height: Math.round(meta.height * scale),
  };
}

async function main() {
  fs.mkdirSync(outdir, { recursive: true });

  const W = 1920;
  const H = 1080;
  const marginX = 52;
  const gap = 28;
  const panelW = Math.floor((W - 2 * marginX - 2 * gap) / 3);
  const panelH = 820;
  const panelY = 176;
  const headerH = 58;
  const captionH = 118;
  const imageH = panelH - headerH - captionH - 34;
  const imageW = panelW - 36;

  const composites = [];
  let svg = `
  <svg width="${W}" height="${H}" xmlns="http://www.w3.org/2000/svg">
    <rect width="${W}" height="${H}" fill="white"/>
    ${textLines(["Neural-network QA: three backup checks people may ask for"], 52, 74, 54, 700, "#0d1424")}
    ${textLines(["The sequence separates ML sanity, phase-space robustness, and fixed-working-point physics risk."], 54, 132, 27, 400, "#565d6a")}
  `;

  for (let i = 0; i < inputs.length; i += 1) {
    const item = inputs[i];
    const x = marginX + i * (panelW + gap);
    const y = panelY;
    svg += `
      <rect x="${x}" y="${y}" width="${panelW}" height="${panelH}" rx="14" fill="#fafbfd" stroke="#dadfe6" stroke-width="2"/>
      ${textLines([item.title], x + 18, y + 43, 30, 700, "#0d1424")}
    `;

    const trimmed = await trimBuffer(item);
    const fit = await fitForPanel(trimmed, imageW, imageH);
    const resized = await sharp(trimmed)
      .resize(fit.width, fit.height, { fit: "inside", withoutEnlargement: true })
      .png()
      .toBuffer();
    const imgX = x + Math.floor((panelW - fit.width) / 2);
    const imgY = y + headerH + 12 + Math.floor((imageH - fit.height) / 2);
    composites.push({ input: resized, left: imgX, top: imgY });

    const capLines = wrapWords(item.caption, 61);
    svg += textLines(capLines, x + 22, y + headerH + imageH + 50, 22, 400, "#0d1424", 1.2);
  }

  svg += textLines(
    [
      "Source: full-stat Au+Au MLP / BDT+MLP validation pack. Stacker plots are backup diagnostics; production promotion still requires runtime parity and closure checks.",
    ],
    54,
    1042,
    18,
    400,
    "#666c77"
  );
  svg += "</svg>";

  const base = await sharp(Buffer.from(svg)).png().toBuffer();
  await sharp(base).composite(composites).png().toFile(outpng);
  console.log(outpng);

  await makeWideAucVersion();
}

async function placePlot(composites, item, x, y, boxW, boxH) {
  const trimmed = await trimBuffer(item);
  const fit = await fitForPanel(trimmed, boxW, boxH);
  const resized = await sharp(trimmed)
    .resize(fit.width, fit.height, { fit: "inside", withoutEnlargement: true })
    .png()
    .toBuffer();
  composites.push({
    input: resized,
    left: x + Math.floor((boxW - fit.width) / 2),
    top: y + Math.floor((boxH - fit.height) / 2),
  });
}

async function makeWideAucVersion() {
  const W = 1920;
  const H = 1080;
  const navy = "#0d1424";
  const gray = "#565d6a";
  const composites = [];

  const top = inputs[1];
  const bottomLeft = inputs[0];
  const bottomRight = inputs[2];

  let svg = `
  <svg width="${W}" height="${H}" xmlns="http://www.w3.org/2000/svg">
    <rect width="${W}" height="${H}" fill="white"/>
    ${textLines(["Neural-network QA backup: does the NN learn, generalize, and operate safely?"], 52, 70, 48, 700, navy)}
    ${textLines(["Three compact checks: training behavior, AUC stability across E_T and centrality, and residual fake rate at WP80."], 54, 122, 25, 400, gray)}

    <rect x="60" y="162" width="1800" height="432" rx="14" fill="#fafbfd" stroke="#dadfe6" stroke-width="2"/>
    ${textLines(["1. Phase-space performance: AUC by E_T and centrality"], 84, 204, 28, 700, navy)}
    ${textLines(["This is the first expert question: does the NN separate signal/background everywhere, or only after integrating bins?"], 84, 236, 19, 400, gray)}

    <rect x="60" y="622" width="875" height="370" rx="14" fill="#fafbfd" stroke="#dadfe6" stroke-width="2"/>
    ${textLines(["2. Training stability"], 84, 662, 26, 700, navy)}
    ${textLines(["Validation flattening marks where extra epochs stop helping unseen candidates."], 84, 690, 18, 400, gray)}

    <rect x="985" y="622" width="875" height="370" rx="14" fill="#fafbfd" stroke="#dadfe6" stroke-width="2"/>
    ${textLines(["3. WP80 operating behavior"], 1009, 662, 26, 700, navy)}
    ${textLines(["The remaining background fake rate identifies the bins that drive purity risk."], 1009, 690, 18, 400, gray)}

    ${textLines(["Source: full-stat Au+Au MLP / BDT+MLP validation pack. Stack plots are backup diagnostics; production promotion still requires runtime parity and closure checks."], 54, 1045, 17, 400, "#666c77")}
  </svg>`;

  await placePlot(composites, top, 92, 248, 1736, 330);
  await placePlot(composites, bottomLeft, 84, 715, 830, 260);
  await placePlot(composites, bottomRight, 1008, 715, 830, 260);

  const base = await sharp(Buffer.from(svg)).png().toBuffer();
  await sharp(base).composite(composites).png().toFile(outpngV2);
  console.log(outpngV2);
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
