import fs from "node:fs/promises";
import fsSync from "node:fs";
import os from "node:os";
import path from "node:path";
import { execFileSync } from "node:child_process";

import {
  Presentation,
  PresentationFile,
  image,
  layers,
  shape,
  text,
} from "@oai/artifact-tool";

const SLIDE_W = 1280;
const SLIDE_H = 720;
const RENDER_SCALE = 2;
const FONT = "Times New Roman";

const C = {
  navy: "#17365D",
  blue: "#2F75B5",
  paleBlue: "#EAF2F8",
  cyan: "#DDEBF7",
  green: "#548235",
  paleGreen: "#E2F0D9",
  red: "#C00000",
  paleRed: "#FCE4D6",
  amber: "#BF7B00",
  paleAmber: "#FFF2CC",
  charcoal: "#263238",
  gray: "#6B7280",
  paleGray: "#F4F6F8",
  border: "#B8C4D2",
  white: "#FFFFFF",
};

function rect(name, x, y, w, h, fill, stroke = fill, radius = 0) {
  return shape({
    name,
    geometry: radius > 0 ? "roundRect" : "rect",
    fill,
    stroke: { color: stroke, width: 1.3 },
    position: { left: x, top: y },
    width: w,
    height: h,
  });
}

function txt(name, value, x, y, w, h, {
  size = 20,
  color = C.charcoal,
  bold = false,
  italic = false,
  align = "left",
  valign = "middle",
  inset = 0,
} = {}) {
  return text([value], {
    name,
    position: { left: x, top: y },
    width: w,
    height: h,
    style: {
      fontSize: `${size}px`,
      typeface: FONT,
      color,
      bold,
      italic,
      alignment: align,
      verticalAlignment: valign,
      autoFit: "shrinkText",
      wrap: "square",
      insets: { top: inset, right: inset, bottom: inset, left: inset },
    },
  });
}

function nodeBox(name, x, y, w, h, extra = {}) {
  return { name, kind: "panel", bbox: [2 * x, 2 * y, 2 * (x + w), 2 * (y + h)], ...extra };
}

function nodeText(name, value, x, y, w, h, size, role = "audience", extra = {}) {
  return {
    name,
    kind: "text",
    role,
    text: value,
    bbox: [2 * x, 2 * y, 2 * (x + w), 2 * (y + h)],
    font_px: 2 * size,
    ...extra,
  };
}

async function saveBlob(blob, outputPath) {
  const bytes = new Uint8Array(await blob.arrayBuffer());
  await fs.writeFile(outputPath, bytes);
}

function eqImage(name, blob, x, y, w, h, alt) {
  return image({
    name,
    blob,
    contentType: "image/png",
    alt,
    fit: "contain",
    position: { left: x, top: y },
    width: w,
    height: h,
  });
}

async function renderLatexAsset(tempDir, name, math, {
  color = C.charcoal,
  paperWidthIn = 8,
  paperHeightIn = 1,
  fontSizePt = 22,
} = {}) {
  const texPath = path.join(tempDir, `${name}.tex`);
  const pdfPath = path.join(tempDir, `${name}.pdf`);
  const pngPrefix = path.join(tempDir, name);
  const pngPath = `${pngPrefix}.png`;
  const lineHeightPt = Math.ceil(fontSizePt * 1.15);
  const tex = String.raw`\documentclass{article}
\usepackage{amsmath}
\usepackage{xcolor}
\usepackage{lmodern}
\usepackage[paperwidth=${paperWidthIn}in,paperheight=${paperHeightIn}in,margin=0in]{geometry}
\pagestyle{empty}
\setlength{\parindent}{0pt}
\begin{document}
\begin{minipage}[c][\paperheight][c]{\paperwidth}
\centering
{\fontsize{${fontSizePt}pt}{${lineHeightPt}pt}\selectfont\color[HTML]{${color.slice(1)}}$\displaystyle ${math}$}
\end{minipage}
\end{document}
`;
  await fs.writeFile(texPath, tex, "utf8");
  execFileSync("pdflatex", [
    "-interaction=nonstopmode",
    "-halt-on-error",
    `-output-directory=${tempDir}`,
    texPath,
  ], { stdio: "pipe" });
  const bundledPdftocairo = "/Users/patsfan753/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/bin/pdftocairo";
  const pdftocairo = fsSync.existsSync(bundledPdftocairo) ? bundledPdftocairo : "pdftocairo";
  execFileSync(pdftocairo, [
    "-png",
    "-transp",
    "-singlefile",
    "-r", "300",
    pdfPath,
    pngPrefix,
  ], { stdio: "pipe" });
  const bundledPython = "/Users/patsfan753/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3";
  const python = fsSync.existsSync(bundledPython) ? bundledPython : "python3";
  execFileSync(python, [
    "-c",
    [
      "from PIL import Image",
      "import sys",
      "p=sys.argv[1]",
      "im=Image.open(p).convert('RGBA')",
      "bbox=im.getchannel('A').getbbox()",
      "im=im.crop(bbox) if bbox else im",
      "pad=8",
      "out=Image.new('RGBA',(im.width+2*pad,im.height+2*pad),(0,0,0,0))",
      "out.paste(im,(pad,pad),im)",
      "out.save(p)",
    ].join(";"),
    pngPath,
  ], { stdio: "pipe" });
  return fs.readFile(pngPath);
}

async function renderMathAssets() {
  const tempDir = await fs.mkdtemp(path.join(os.tmpdir(), "purity-slide-math-"));
  return {
    abcd: await renderLatexAsset(
      tempDir,
      "abcd",
      String.raw`S_A = A - \frac{(B-f_B S_A)(C-f_C S_A)}{D-f_D S_A}`,
      { paperWidthIn: 10.4, paperHeightIn: 1.0, fontSizePt: 24 },
    ),
    purity: await renderLatexAsset(
      tempDir,
      "purity",
      String.raw`\color[HTML]{548235} P_{\gamma}(p_T^{\gamma})=\frac{S_A}{A}\qquad \color[HTML]{C00000} N_{\mathrm{bkg}}^{A}=A-S_A`,
      { paperWidthIn: 10.3, paperHeightIn: 1.0, fontSizePt: 24 },
    ),
    sideband: await renderLatexAsset(
      tempDir,
      "sideband",
      String.raw`\begin{aligned}N_{\mathrm{bkg}}^{C} &= C-f_C S_A\\[2pt]\alpha_C &= \frac{N_{\mathrm{bkg}}^{A}}{N_{\mathrm{bkg}}^{C}}\end{aligned}`,
      { paperWidthIn: 4.4, paperHeightIn: 1.5, fontSizePt: 21 },
    ),
    correction: await renderLatexAsset(
      tempDir,
      "correction",
      String.raw`H_A^{\mathrm{sig}}(x_{J\gamma})=\frac{H_A(x_{J\gamma})-\alpha_C H_C(x_{J\gamma})}{1-\alpha_C f_C}`,
      { color: C.navy, paperWidthIn: 6.8, paperHeightIn: 1.2, fontSizePt: 20 },
    ),
    hist: await renderLatexAsset(
      tempDir,
      "hist",
      String.raw`H_C(x_{J\gamma})`,
      { color: C.blue, paperWidthIn: 3.4, paperHeightIn: 1.0, fontSizePt: 20 },
    ),
    xj: await renderLatexAsset(
      tempDir,
      "xj",
      String.raw`x_{J\gamma}`,
      { paperWidthIn: 2.3, paperHeightIn: 1.0, fontSizePt: 15 },
    ),
    final: await renderLatexAsset(
      tempDir,
      "final",
      String.raw`\frac{1}{N_{\gamma}}\frac{dN_{\mathrm{jet}}}{dx_{J\gamma}}=\frac{H_{A,\mathrm{unfolded}}^{\mathrm{sig}}(x_{J\gamma})}{S_{A}^{\mathrm{unfolded}}\,\Delta x_{J\gamma}}`,
      { color: C.navy, paperWidthIn: 13.2, paperHeightIn: 1.0, fontSizePt: 23 },
    ),
  };
}

function buildSlide(presentation, mathAssets) {
  const slide = presentation.slides.add();
  const elements = [];
  const nodes = [];

  elements.push(rect("canvas", 0, 0, SLIDE_W, SLIDE_H, C.white, C.white));

  elements.push(txt(
    "title",
    "Purity corrects the photon tag — not the recoil requirement",
    44, 24, 1120, 54,
    { size: 39, color: C.navy, bold: true, valign: "top" },
  ));
  nodes.push(nodeText(
    "title",
    "Purity corrects the photon tag — not the recoil requirement",
    44, 24, 1120, 54, 39, "title",
  ));

  elements.push(txt(
    "subtitle",
    "Use the same event-leading photon population as the p+p recoil measurement, in each photon transverse-momentum bin.",
    46, 92, 1148, 30,
    { size: 21, color: C.gray, valign: "top" },
  ));
  nodes.push(nodeText(
    "subtitle",
    "Use the same event-leading photon population as the p+p recoil measurement, in each photon transverse-momentum bin.",
    46, 92, 1148, 30, 21,
  ));

  const left = { x: 44, y: 128, w: 325, h: 456 };
  const mid = { x: 392, y: 128, w: 484, h: 456 };
  const right = { x: 899, y: 128, w: 337, h: 456 };
  for (const [name, box] of [["left panel", left], ["middle panel", mid], ["right panel", right]]) {
    elements.push(rect(name, box.x, box.y, box.w, box.h, C.white, C.border, 10));
    nodes.push(nodeBox(name, box.x, box.y, box.w, box.h));
  }

  elements.push(rect("left step band", 44, 128, 325, 52, C.paleBlue, C.paleBlue, 10));
  elements.push(txt("left step", "1   Count photon tags with ABCD", 60, 136, 292, 36,
    { size: 23, color: C.navy, bold: true }));
  nodes.push(nodeText("left step", "1   Count photon tags with ABCD", 60, 136, 292, 36, 23));

  const gridX = 150;
  const gridY = 242;
  const cellW = 94;
  const cellH = 88;
  const gap = 8;
  elements.push(txt("tight column", "Tight", gridX, 204, cellW, 28,
    { size: 20, color: C.charcoal, bold: true, align: "center" }));
  elements.push(txt("non-tight column", "Non-tight", gridX + cellW + gap, 204, cellW, 28,
    { size: 20, color: C.charcoal, bold: true, align: "center" }));
  elements.push(txt("isolated row", "Isolated", 50, gridY + 22, 94, 45,
    { size: 18, color: C.charcoal, bold: true, align: "center" }));
  elements.push(txt("nonisolated row", "Non-isolated", 50, gridY + cellH + gap + 9, 94, 64,
    { size: 17, color: C.charcoal, bold: true, align: "center" }));
  nodes.push(nodeText("tight column", "Tight", gridX, 204, cellW, 28, 20));
  nodes.push(nodeText("non-tight column", "Non-tight", gridX + cellW + gap, 204, cellW, 28, 20));
  nodes.push(nodeText("isolated row", "Isolated", 50, gridY + 22, 94, 45, 18));
  nodes.push(nodeText("nonisolated row", "Non-isolated", 50, gridY + cellH + gap + 9, 94, 64, 17));

  const cells = [
    ["A cell", "A", gridX, gridY, C.paleGreen, C.green, "signal-tag region"],
    ["C cell", "C", gridX + cellW + gap, gridY, C.cyan, C.blue, "fails at least 2 tight cuts"],
    ["B cell", "B", gridX, gridY + cellH + gap, C.paleGray, C.gray, "ABCD control"],
    ["D cell", "D", gridX + cellW + gap, gridY + cellH + gap, C.paleGray, C.gray, "ABCD control"],
  ];
  for (const [name, label, x, y, fill, stroke, sub] of cells) {
    elements.push(rect(name, x, y, cellW, cellH, fill, stroke, 8));
    elements.push(txt(`${name} label`, label, x + 8, y + 7, cellW - 16, 38,
      { size: 31, color: stroke, bold: true, align: "center" }));
    elements.push(txt(`${name} sub`, sub, x + 5, y + 47, cellW - 10, 27,
      { size: 14, color: C.charcoal, align: "center" }));
    nodes.push(nodeText(`${name} label`, label, x + 8, y + 7, cellW - 16, 38, 31));
    nodes.push(nodeText(`${name} sub`, sub, x + 5, y + 47, cellW - 10, 27, 14));
  }

  elements.push(txt(
    "tag population",
    "A, B, C, D each count one event-leading photon tag per event and region.",
    64, 448, 286, 58,
    { size: 18, color: C.charcoal, bold: true, align: "center" },
  ));
  elements.push(txt(
    "isolation gap note",
    "The isolation gap is excluded from the non-isolated sideband.",
    68, 507, 278, 34,
    { size: 16, color: C.gray, italic: true, align: "center" },
  ));
  elements.push(txt(
    "no jet count",
    "No recoil jet is required for the purity counts.",
    72, 542, 270, 32,
    { size: 19, color: C.red, bold: true, align: "center" },
  ));
  nodes.push(nodeText(
    "tag population",
    "A, B, C, D each count one event-leading photon tag per event and region.",
    64, 448, 286, 58, 18,
  ));
  nodes.push(nodeText(
    "isolation gap note",
    "The isolation gap is excluded from the non-isolated sideband.",
    68, 507, 278, 34, 16,
  ));
  nodes.push(nodeText(
    "no jet count",
    "No recoil jet is required for the purity counts.",
    72, 542, 270, 32, 19,
  ));

  elements.push(rect("middle step band", 392, 128, 484, 52, C.paleBlue, C.paleBlue, 10));
  elements.push(txt("middle step", "2   Solve signal tags, then subtract Region C", 408, 136, 450, 36,
    { size: 23, color: C.navy, bold: true }));
  nodes.push(nodeText("middle step", "2   Solve signal tags, then subtract Region C", 408, 136, 450, 36, 23));

  elements.push(rect("abcd equation band", 412, 194, 444, 84, C.paleAmber, C.amber, 7));
  elements.push(eqImage(
    "abcd equation",
    mathAssets.abcd,
    426, 200, 416, 50,
    "Leakage-corrected ABCD signal-count equation",
  ));
  elements.push(txt(
    "abcd implementation note",
    "Leakage corrected in B, C, and D; current background factor = 1",
    430, 247, 408, 22,
    { size: 17, color: C.gray, italic: true, align: "center" },
  ));
  nodes.push(nodeText(
    "abcd equation",
    "S_A = A − [(B − f_B S_A)(C − f_C S_A)] / (D − f_D S_A)",
    426, 207, 416, 40, 23,
  ));
  nodes.push(nodeText(
    "abcd implementation note",
    "Leakage corrected in B, C, and D; current background factor = 1",
    430, 247, 408, 22, 17,
  ));

  elements.push(eqImage(
    "purity and background counts",
    mathAssets.purity,
    420, 279, 432, 50,
    "Photon purity and Region-A background-count equations",
  ));
  nodes.push(nodeText("purity equation", "Pγ(pTγ) = S_A / A", 420, 286, 215, 36, 25));
  nodes.push(nodeText("background count", "N_bkg^A = A − S_A", 638, 286, 214, 36, 23));

  elements.push(txt("region C heading", "Region C provides the fake-tag recoil shape", 414, 329, 438, 30,
    { size: 21, color: C.blue, bold: true, align: "center" }));
  nodes.push(nodeText("region C heading", "Region C provides the fake-tag recoil shape", 414, 329, 438, 30, 21));

  const histX = 430;
  const histY = 387;
  const histW = 126;
  const histH = 129;
  elements.push(rect("hist background", histX, histY, histW, histH, C.white, C.border, 3));
  elements.push(rect("hist x axis", histX + 14, histY + histH - 20, histW - 25, 2, C.charcoal, C.charcoal));
  elements.push(rect("hist y axis", histX + 14, histY + 12, 2, histH - 31, C.charcoal, C.charcoal));
  const barHeights = [34, 67, 86, 74, 51, 28];
  for (let i = 0; i < barHeights.length; i += 1) {
    const bw = 11;
    const bx = histX + 22 + i * 15;
    const bh = barHeights[i];
    elements.push(rect(`hist bar ${i + 1}`, bx, histY + histH - 21 - bh, bw, bh, C.cyan, C.blue));
  }
  elements.push(eqImage(
    "hist label",
    mathAssets.hist,
    histX + 20, histY + 5, histW - 30, 30,
    "Region-C recoil histogram label",
  ));
  elements.push(eqImage(
    "hist x label",
    mathAssets.xj,
    histX + 70, histY + histH - 21, 45, 21,
    "x J gamma axis label",
  ));
  elements.push(txt("schematic label", "schematic", histX + 23, histY + histH + 2, histW - 34, 18,
    { size: 14, color: C.gray, italic: true, align: "center" }));
  nodes.push(nodeText("hist label", "H_C(xJ)", histX + 20, histY + 8, histW - 30, 24, 19, "plot_annotation"));
  nodes.push(nodeText("hist x label", "xJγ", histX + 70, histY + histH - 19, 45, 17, 15, "plot_annotation"));
  nodes.push(nodeText("schematic label", "schematic", histX + 23, histY + histH + 2, histW - 34, 18, 14, "plot_annotation"));

  elements.push(eqImage(
    "sideband normalization equations",
    mathAssets.sideband,
    582, 372, 260, 78,
    "Region-C background count and normalization equations",
  ));
  elements.push(rect("hist correction band", 576, 456, 272, 75, C.paleBlue, C.blue, 7));
  elements.push(eqImage(
    "hist correction equation",
    mathAssets.correction,
    588, 461, 248, 58,
    "Leakage-corrected Region-C recoil subtraction equation",
  ));
  elements.push(txt(
    "shape normalization note",
    "Region C gives the shape; the ABCD normalization fixes how much to subtract.",
    570, 539, 282, 36,
    { size: 17, color: C.gray, italic: true, align: "center" },
  ));
  nodes.push(nodeText("c background", "N_bkg^C = C − f_C S_A", 582, 378, 260, 32, 22));
  nodes.push(nodeText("alpha equation", "α_C = N_bkg^A / N_bkg^C", 582, 415, 260, 32, 22));
  nodes.push(nodeText(
    "hist correction equation",
    "H_A^sig(xJ) = [H_A(xJ) − α_C H_C(xJ)] / [1 − α_C f_C]",
    588, 466, 248, 48, 20,
  ));
  nodes.push(nodeText(
    "shape normalization note",
    "Region C gives the shape; the ABCD normalization fixes how much to subtract.",
    570, 539, 282, 36, 17,
  ));

  elements.push(rect("right step band", 899, 128, 337, 52, C.paleBlue, C.paleBlue, 10));
  elements.push(txt("right step", "3   Match the analysis population", 915, 136, 305, 36,
    { size: 23, color: C.navy, bold: true }));
  nodes.push(nodeText("right step", "3   Match the analysis population", 915, 136, 305, 36, 23));

  elements.push(rect("use band", 921, 198, 293, 148, C.paleGreen, C.green, 8));
  elements.push(txt("use label", "USE", 937, 208, 62, 28,
    { size: 22, color: C.green, bold: true, align: "center" }));
  elements.push(txt(
    "use text",
    "• One event-leading photon tag per event and ABCD region\n• Same photon-momentum bin and cuts as the recoil measurement\n• No away-side recoil requirement",
    938, 240, 260, 94,
    { size: 19, color: C.charcoal, valign: "top" },
  ));
  nodes.push(nodeText("use label", "USE", 937, 208, 62, 28, 22));
  nodes.push(nodeText(
    "use text",
    "• One event-leading photon tag per event and ABCD region\n• Same photon-momentum bin and cuts as the recoil measurement\n• No away-side recoil requirement",
    938, 240, 260, 94, 19,
  ));

  elements.push(rect("region C contract band", 921, 360, 293, 91, C.cyan, C.blue, 8));
  elements.push(txt("region C contract label", "FOR THE REGION-C HISTOGRAM", 937, 369, 245, 26,
    { size: 20, color: C.blue, bold: true }));
  elements.push(txt(
    "region C contract text",
    "Use one event-leading isolated ∧ non-tight C tag, then pair it with all selected away-side recoil jets — not only the leading jet.",
    938, 397, 258, 45,
    { size: 16, color: C.charcoal, valign: "top" },
  ));
  nodes.push(nodeText("region C contract label", "FOR THE REGION-C HISTOGRAM", 937, 369, 245, 26, 20));
  nodes.push(nodeText(
    "region C contract text",
    "Use one event-leading isolated ∧ non-tight C tag, then pair it with all selected away-side recoil jets — not only the leading jet.",
    938, 397, 258, 45, 16,
  ));

  elements.push(rect("do not band", 921, 465, 293, 99, C.paleRed, C.red, 8));
  elements.push(txt("do not label", "DO NOT", 937, 475, 91, 27,
    { size: 21, color: C.red, bold: true, align: "center" }));
  elements.push(txt(
    "do not text",
    "• Condition photon purity on finding a recoil jet\n• Substitute all-candidate purity without a leading-tag transfer check",
    938, 505, 258, 50,
    { size: 17, color: C.charcoal, valign: "top" },
  ));
  nodes.push(nodeText("do not label", "DO NOT", 937, 475, 91, 27, 21));
  nodes.push(nodeText(
    "do not text",
    "• Condition photon purity on finding a recoil jet\n• Substitute all-candidate purity without a leading-tag transfer check",
    938, 505, 258, 50, 17,
  ));

  elements.push(rect("bottom conclusion band", 44, 603, 1192, 87, C.paleBlue, C.blue, 9));
  elements.push(eqImage(
    "bottom equation",
    mathAssets.final,
    120, 600, 1040, 58,
    "Final unfolded per-photon differential recoil-yield equation",
  ));
  elements.push(txt(
    "bottom takeaway",
    "Every accepted leading photon belongs in the denominator — even when no recoil jet is found.",
    88, 658, 1104, 24,
    { size: 22, color: C.red, bold: true, align: "center" },
  ));
  nodes.push(nodeText(
    "bottom equation",
    "(1/Nγ) dNjet/dxJ  =  H_A^sig,unfolded(xJ) / [S_A^unfolded · ΔxJ]",
    72, 611, 1136, 36, 27,
  ));
  nodes.push(nodeText(
    "bottom takeaway",
    "Every accepted leading photon belongs in the denominator — even when no recoil jet is found.",
    88, 658, 1104, 24, 22,
  ));

  slide.compose(layers({ name: "purity-correction-contract", width: "fill", height: "fill" }, elements));
  return { slide, nodes };
}

async function main() {
  const outputDir = path.resolve(process.argv[2] || ".");
  await fs.mkdir(outputDir, { recursive: true });

  const stem = "slide_purity_correction_regionC_contract_v1";
  const pptxPath = path.join(outputDir, `${stem}.pptx`);
  const pngPath = path.join(outputDir, `${stem}.png`);
  const layoutPath = path.join(outputDir, `${stem}_layout_nodes.json`);
  const artifactLayoutPath = path.join(outputDir, `${stem}_artifact_layout.json`);
  const manifestPath = path.join(outputDir, `${stem}_manifest.json`);
  const speakerScriptPath = path.join(outputDir, `${stem}_speaker_script.md`);

  const mathAssets = await renderMathAssets();
  const presentation = Presentation.create({ slideSize: { width: SLIDE_W, height: SLIDE_H } });
  const { slide, nodes } = buildSlide(presentation, mathAssets);

  const preview = await presentation.export({ slide, format: "png", scale: RENDER_SCALE });
  await saveBlob(preview, pngPath);
  const bundledPython = "/Users/patsfan753/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3";
  const python = fsSync.existsSync(bundledPython) ? bundledPython : "python3";
  execFileSync(python, [
    "-c",
    "from PIL import Image; import sys; p=sys.argv[1]; Image.open(p).convert('RGB').save(p)",
    pngPath,
  ], { stdio: "pipe" });
  const artifactLayout = await presentation.export({ slide, format: "layout" });
  await saveBlob(artifactLayout, artifactLayoutPath);
  const pptx = await PresentationFile.exportPptx(presentation);
  await pptx.save(pptxPath);

  const layoutPayload = {
    slide_size: [SLIDE_W * RENDER_SCALE, SLIDE_H * RENDER_SCALE],
    title_axis_x: 88,
    minimum_audience_font_px: 28,
    minimum_title_font_px: 52,
    minimum_plot_annotation_font_px: 28,
    nodes,
  };
  await fs.writeFile(layoutPath, `${JSON.stringify(layoutPayload, null, 2)}\n`, "utf8");

  const manifest = {
    artifact_id: "the87_purity_correction_regionC_contract_v1",
    title: "Purity corrects the photon tag — not the recoil requirement",
    canvas: "16:9, 2560x1440 PNG",
    outputs: {
      png: pngPath,
      pptx: pptxPath,
      layout_nodes: layoutPath,
      artifact_layout: artifactLayoutPath,
      speaker_script: speakerScriptPath,
    },
    science_contract: {
      purity_population: "One event-leading photon candidate per event and ABCD region, in the same pTgamma bin and photon selection used by xJgamma, without requiring an away-side recoil jet.",
      region_c_histogram: "H_C(xJ) uses one event-leading isolated non-tight Region C tag, failing at least two tight cuts, paired with all selected away-side recoil jets; Region C supplies the fake-tag recoil shape.",
      current_abcd_equation: "S_A = A - (B-f_B*S_A)(C-f_C*S_A)/(D-f_D*S_A), with R_bkg=1 in the current implementation.",
      purity: "Pgamma = S_A/A and N_bkg^A = A-S_A.",
      region_c_correction: "N_bkg^C=C-f_C*S_A, alpha_C=N_bkg^A/N_bkg^C, H_A^sig=(H_A-alpha_C*H_C)/(1-alpha_C*f_C).",
      final_normalization: "Unfold the corrected H_A^sig and divide by the unfolded photon signal count S_A and the xJ bin width Delta xJ.",
    },
    source_evidence: [
      "src/RecoilJets.cc event-leading xJ purity counters and Region C recoil filling",
      "macros/AnalyzeRecoilJets_RooUnfoldPipeline.cpp ComputeABCDSignalCounts and leakage-corrected Region C subtraction",
      "macros/AnalyzeRecoilJets.h SolveLeakageCorrectedSA",
    ],
    caveats: [
      "This slide states the current xJgamma tag-level contract; candidate-level purity is not interchangeable unless transfer to the event-leading tag population is demonstrated.",
      "The histogram is schematic and contains no data points.",
      "Google Slides was not mutated.",
    ],
    speaker_script: [
      "Purity is the fraction of signal photons in the same leading-tag population that defines the measurement, not the fraction after asking for a recoil jet.",
      "The ABCD counts are therefore filled once per event and region before recoil matching; the signal count is S_A and the purity is S_A divided by A.",
      "Region C has a different role in the numerator: after selecting one event-leading isolated non-tight Region C tag, its histogram of all selected away-side recoil jets supplies the fake-tag recoil shape, normalized by alpha_C and corrected for signal leakage.",
      "Finally we unfold that corrected recoil spectrum and divide by the unfolded signal-photon count and the xJ bin width; requiring a recoil jet in the purity denominator would bias the per-photon yield.",
    ],
  };
  await fs.writeFile(manifestPath, `${JSON.stringify(manifest, null, 2)}\n`, "utf8");
  await fs.writeFile(
    speakerScriptPath,
    [
      "# Purity correction and Region C — speaker script",
      "",
      "Purity is the fraction of signal photons in the same event-leading tag population that defines the xJgamma measurement. It is not the fraction after asking for a recoil jet.",
      "",
      "The ABCD counts are filled once per event and region before recoil matching. The leakage-corrected solution gives the signal-tag count S_A, and the photon purity is S_A divided by A.",
      "",
      "Region C has a different role in the recoil numerator. After selecting the event-leading isolated non-tight tag, which fails at least two tight cuts, we pair it with all selected away-side recoil jets to build H_C(xJ). That histogram provides the fake-tag recoil shape; alpha_C supplies its normalization, and the denominator in the correction removes prompt-photon leakage into Region C. The isolation gap is excluded from the non-isolated sideband.",
      "",
      "Finally, we unfold the corrected recoil spectrum and divide by the unfolded signal-photon count and the xJ bin width. Every accepted leading photon belongs in that denominator, including events where no recoil jet is found. Conditioning the purity on a recoil jet would bias the per-photon yield.",
      "",
    ].join("\n"),
    "utf8",
  );

  console.log(JSON.stringify({
    pptxPath,
    pngPath,
    layoutPath,
    artifactLayoutPath,
    manifestPath,
    speakerScriptPath,
  }, null, 2));
}

main().catch((error) => {
  console.error(error.stack || String(error));
  process.exit(1);
});
