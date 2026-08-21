import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { execFileSync } from "node:child_process";

const artifactToolModule = process.env.ARTIFACT_TOOL_MODULE || "@oai/artifact-tool";
const {
  Presentation,
  PresentationFile,
  image,
  layers,
  shape,
  text,
} = await import(artifactToolModule);

const SLIDE_W = 1280;
const SLIDE_H = 720;
const SCALE = 2;
const FONT = "Times New Roman";

const C = {
  navy: "#17365D",
  blue: "#2F75B5",
  teal: "#2F7F85",
  paleBlue: "#EAF2F8",
  paleTeal: "#E8F3F3",
  green: "#548235",
  paleGreen: "#E2F0D9",
  amber: "#B7791F",
  paleAmber: "#FFF4D6",
  red: "#B91C1C",
  charcoal: "#202124",
  gray: "#5F6670",
  line: "#B8C4D2",
  panel: "#FBFCFD",
  stage: "#F2F5F8",
  white: "#FFFFFF",
};

function rect(name, x, y, w, h, fill, stroke = fill, radius = 8, width = 1.2) {
  return shape({
    name,
    geometry: radius > 0 ? "roundRect" : "rect",
    fill,
    stroke: { color: stroke, width },
    position: { left: x, top: y },
    width: w,
    height: h,
  });
}

function arrow(name, x, y, w, h, direction = "right", color = C.blue) {
  return shape({
    name,
    geometry: direction === "down" ? "downArrow" : "rightArrow",
    fill: color,
    stroke: { color, width: 0.8 },
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
  inset = 8,
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

function eq(name, blob, x, y, w, h, alt) {
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

function layoutNode(name, kind, x, y, w, h, extra = {}) {
  return {
    name,
    kind,
    bbox: [SCALE * x, SCALE * y, SCALE * (x + w), SCALE * (y + h)],
    ...extra,
  };
}

async function saveBlob(blob, outputPath) {
  await fs.writeFile(outputPath, new Uint8Array(await blob.arrayBuffer()));
}

async function renderLatexAsset(tempDir, name, math, {
  widthIn = 8,
  heightIn = 0.8,
  fontSizePt = 22,
  color = C.charcoal,
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
\usepackage[paperwidth=${widthIn}in,paperheight=${heightIn}in,margin=0in]{geometry}
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
  const pdftocairoBin = process.env.PDFTOCAIRO_BIN || "pdftocairo";
  execFileSync(pdftocairoBin, [
    "-png", "-transp", "-singlefile", "-r", "300", pdfPath, pngPrefix,
  ], { stdio: "pipe" });
  const pythonBin = process.env.PYTHON_BIN || "python3";
  execFileSync(pythonBin, [
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

async function renderMath() {
  const tempDir = await fs.mkdtemp(path.join(os.tmpdir(), "the219-purity-flow-math-"));
  return {
    leakage: await renderLatexAsset(
      tempDir,
      "leakage",
      String.raw`S_A=A-\frac{(B-f_BS_A)(C-f_CS_A)}{D-f_DS_A}`,
      { widthIn: 8.6, heightIn: 0.95, fontSizePt: 29, color: C.navy },
    ),
    signalCount: await renderLatexAsset(
      tempDir,
      "signal_count",
      String.raw`S_A`,
      { widthIn: 1.4, heightIn: 0.65, fontSizePt: 29, color: C.navy },
    ),
    fakeCount: await renderLatexAsset(
      tempDir,
      "fake_count",
      String.raw`N_{\mathrm{bkg}}^{A}=A-S_A`,
      { widthIn: 4.1, heightIn: 0.7, fontSizePt: 27, color: C.navy },
    ),
    regionA: await renderLatexAsset(
      tempDir,
      "regionA",
      String.raw`H_A(x_{J\gamma})`,
      { widthIn: 3.2, heightIn: 0.75, fontSizePt: 34, color: C.navy },
    ),
    fake: await renderLatexAsset(
      tempDir,
      "fake",
      String.raw`H_{\mathrm{fake}}^{A}(x_{J\gamma})=N_{\mathrm{bkg}}^{A}\left[\frac{H_C(x_{J\gamma})}{N_C}\right]`,
      { widthIn: 7.1, heightIn: 1.0, fontSizePt: 29, color: C.navy },
    ),
    final: await renderLatexAsset(
      tempDir,
      "final",
      String.raw`H_{\mathrm{signal}}(x_{J\gamma})=H_A(x_{J\gamma})-H_{\mathrm{fake}}^{A}(x_{J\gamma})`,
      { widthIn: 8.8, heightIn: 0.9, fontSizePt: 32, color: C.navy },
    ),
  };
}

function buildSlide(presentation, math) {
  const slide = presentation.slides.add();
  const elements = [];
  const nodes = [];

  elements.push(rect("canvas", 0, 0, SLIDE_W, SLIDE_H, C.white, C.white, 0, 0));

  // Connectors are composed first so every arrow stays behind the content.

  elements.push(txt(
    "title",
    "Purity Correction Method Choice",
    34, 16, 1185, 55,
    { size: 40, color: "#000000", bold: true, valign: "top", inset: 0 },
  ));
  nodes.push(layoutNode("title", "text", 34, 16, 1185, 55, {
    role: "title", text: "Purity Correction Method Choice", font_px: 80,
  }));
  elements.push(txt(
    "method label",
    "ATLAS: two-stage subtraction",
    930, 35, 302, 28,
    { size: 21, color: C.navy, bold: true, align: "right", inset: 0 },
  ));

  const stage1 = { x: 46, y: 82, w: 1188, h: 245 };
  elements.push(rect("stage1 field", stage1.x, stage1.y, stage1.w, stage1.h, "#EDF5FB", "#AFCDE5", 12, 1.2));
  elements.push(txt(
    "stage1 question",
    "1. How many fake photons are in Region A?",
    65, 92, 560, 34,
    { size: 28, color: C.navy, bold: true, inset: 0 },
  ));
  nodes.push(layoutNode("stage1 field", "panel", stage1.x, stage1.y, stage1.w, stage1.h));

  const abcd = { x: 65, y: 130, w: 348, h: 177 };
  const correction = { x: 450, y: 130, w: 765, h: 177 };
  elements.push(rect("ABCD map", abcd.x, abcd.y, abcd.w, abcd.h, C.white, "#87B4D8", 9, 1.2));
  elements.push(rect("Correction calculation", correction.x, correction.y, correction.w, correction.h, C.white, "#87B4D8", 9, 1.2));
  nodes.push(layoutNode("ABCD map", "panel", abcd.x, abcd.y, abcd.w, abcd.h));
  nodes.push(layoutNode("Correction calculation", "panel", correction.x, correction.y, correction.w, correction.h));

  elements.push(txt("ABCD title", "Standard ABCD counts", 80, 137, 318, 27,
    { size: 25, color: C.navy, bold: true, align: "center", inset: 0 }));
  elements.push(txt("isolated label", "isolated", 165, 166, 104, 20,
    { size: 18, color: C.gray, bold: true, align: "center", inset: 0 }));
  elements.push(txt("nonisolated label", "non-isolated", 274, 166, 111, 20,
    { size: 18, color: C.gray, bold: true, align: "center", inset: 0 }));
  elements.push(txt("tight label", "tight", 81, 199, 68, 34,
    { size: 19, color: C.gray, bold: true, align: "right", inset: 0 }));
  elements.push(txt("nontight label", "non-tight", 73, 242, 76, 34,
    { size: 19, color: C.gray, bold: true, align: "right", inset: 0 }));
  const cells = [
    ["A cell", 158, 193, "#E2F0D9", C.green, "A  •  signal"],
    ["B cell", 272, 193, "#F2F5F8", C.line, "B"],
    ["C cell", 158, 238, "#F2F5F8", C.line, "C"],
    ["D cell", 272, 238, "#F2F5F8", C.line, "D"],
  ];
  for (const [name, x, y, fill, stroke, label] of cells) {
    elements.push(rect(name, x, y, 108, 40, fill, stroke, 5, 1.1));
    elements.push(txt(`${name} label`, label, x, y, 108, 40,
      { size: 21, color: C.navy, bold: true, align: "center", inset: 0 }));
  }
  elements.push(txt("ABCD foot", "A = tight + isolated signal region", 92, 283, 294, 18,
    { size: 18, color: C.charcoal, bold: true, align: "center", inset: 0 }));
  elements.push(txt("ABCD flow arrow", "→", 414, 195, 36, 40,
    { size: 34, color: C.blue, bold: true, align: "center", inset: 0 }));

  elements.push(txt("Correction title", "Leakage-corrected prompt count", 470, 138, 725, 28,
    { size: 26, color: C.navy, bold: true, align: "center", inset: 0 }));
  elements.push(txt(
    "Correction explanation",
    "MC estimates genuine-photon leakage A → B, C, D  (same correction logic as PPG12)",
    480, 169, 705, 25,
    { size: 20, color: C.charcoal, align: "center", inset: 0 },
  ));
  elements.push(eq("Leakage equation", math.leakage, 523, 194, 620, 62,
    "Leakage-corrected ABCD signal-count equation"));

  elements.push(rect("Signal count result", 480, 259, 327, 37, "#E8F3F3", "#9EC6C8", 6, 1));
  elements.push(eq("Signal count symbol", math.signalCount, 495, 264, 58, 27,
    "Genuine prompt-photon count in Region A"));
  elements.push(txt("Signal count meaning", "genuine prompt photons in A", 554, 262, 238, 31,
    { size: 19, color: C.charcoal, bold: true, align: "center", inset: 0 }));
  elements.push(rect("Fake count result", 826, 259, 359, 37, "#FFF0D5", "#D9B368", 6, 1));
  elements.push(eq("Fake count symbol", math.fakeCount, 842, 263, 154, 29,
    "Fake-photon count in Region A"));
  elements.push(txt("Fake count meaning", "fake photons in A", 997, 262, 170, 31,
    { size: 19, color: C.charcoal, bold: true, align: "center", inset: 0 }));

  const stage2 = { x: 46, y: 346, w: 1188, h: 338 };
  elements.push(rect("stage2 field", stage2.x, stage2.y, stage2.w, stage2.h, "#FFF7E8", "#E5C88A", 12, 1.2));
  elements.push(txt(
    "stage2 question",
    "2. Remove recoil jets carried by fake photon tags",
    65, 357, 690, 34,
    { size: 27, color: C.navy, bold: true, inset: 0 },
  ));
  elements.push(txt(
    "stage2 handoff",
    "Shape from Region C  •  normalization from Stage 1",
    743, 360, 472, 28,
    { size: 20, color: C.amber, bold: true, align: "right", inset: 0 },
  ));
  nodes.push(layoutNode("stage2 field", "panel", stage2.x, stage2.y, stage2.w, stage2.h));

  const recoA = { x: 65, y: 410, w: 235, h: 232 };
  const fakeTemplate = { x: 350, y: 410, w: 430, h: 232 };
  const corrected = { x: 835, y: 410, w: 380, h: 232 };
  elements.push(rect("Region A recoil", recoA.x, recoA.y, recoA.w, recoA.h, C.white, "#7EA7CA", 9, 1.3));
  elements.push(rect("Fake recoil template", fakeTemplate.x, fakeTemplate.y, fakeTemplate.w, fakeTemplate.h, C.white, "#D4AA58", 9, 1.3));
  elements.push(rect("Purity corrected recoil", corrected.x, corrected.y, corrected.w, corrected.h, C.white, "#8DB47A", 9, 1.3));
  elements.push(rect("Region A rail", recoA.x, recoA.y, recoA.w, 7, C.blue, C.blue, 0, 0));
  elements.push(rect("Template rail", fakeTemplate.x, fakeTemplate.y, fakeTemplate.w, 7, C.amber, C.amber, 0, 0));
  elements.push(rect("Corrected rail", corrected.x, corrected.y, corrected.w, 7, C.green, C.green, 0, 0));
  nodes.push(layoutNode("Region A recoil", "panel", recoA.x, recoA.y, recoA.w, recoA.h));
  nodes.push(layoutNode("Fake recoil template", "panel", fakeTemplate.x, fakeTemplate.y, fakeTemplate.w, fakeTemplate.h));
  nodes.push(layoutNode("Purity corrected recoil", "panel", corrected.x, corrected.y, corrected.w, corrected.h));

  elements.push(txt("Region A title", "Region A recoil", 80, 428, 205, 30,
    { size: 26, color: C.navy, bold: true, align: "center", inset: 0 }));
  elements.push(eq("Region A equation", math.regionA, 90, 468, 185, 63,
    "Region A xJgamma recoil distribution"));
  elements.push(txt("Region A contents", "genuine photon+jet\n+\nfake-photon/dijet", 84, 540, 197, 78,
    { size: 21, color: C.charcoal, bold: true, align: "center", valign: "top", inset: 0 }));

  elements.push(txt("Template title", "Fake-recoil template from Region C", 366, 428, 398, 30,
    { size: 25, color: C.navy, bold: true, align: "center", inset: 0 }));
  elements.push(txt(
    "Region C definition",
    "isolated + non-tight; dominated by isolated neutral mesons\nand similar jet fragments that resemble photons",
    375, 468, 380, 50,
    { size: 20, color: C.charcoal, align: "center", valign: "top", inset: 0 },
  ));
  elements.push(eq("Scaled fake equation", math.fake, 378, 520, 374, 75,
    "Region C per-candidate recoil shape scaled to the fake-photon count in Region A"));
  elements.push(txt(
    "Template meaning",
    "per-candidate C shape  ×  fake-photon count in A",
    373, 602, 384, 25,
    { size: 19, color: C.charcoal, bold: true, align: "center", inset: 0 },
  ));

  elements.push(txt("minus operator", "−", 304, 488, 42, 70,
    { size: 62, color: C.navy, bold: false, align: "center", inset: 0 }));
  elements.push(txt("equals operator", "=", 790, 488, 42, 70,
    { size: 48, color: C.green, bold: true, align: "center", inset: 0 }));

  elements.push(txt("Corrected title", "Purity-corrected recoil", 851, 428, 348, 30,
    { size: 26, color: C.navy, bold: true, align: "center", inset: 0 }));
  elements.push(eq("Final subtraction", math.final, 855, 474, 340, 88,
    "ATLAS fake-photon recoil subtraction equation"));
  elements.push(txt(
    "Final meaning",
    "subtract the scaled Region-C template\nfrom the measured Region-A distribution",
    866, 575, 318, 54,
    { size: 21, color: C.charcoal, bold: true, align: "center", valign: "top", inset: 0 },
  ));

  slide.compose(layers({ name: "the219-purity-method-flow", width: "fill", height: "fill" }, elements));
  slide.speakerNotes.textFrame.setText([
    "ATLAS treats the purity correction as two linked subtractions.",
    "First, the standard tight or non-tight and isolated or non-isolated ABCD regions determine how many selected Region-A candidates are genuine photons. Simulation supplies the signal leakage fractions into B, C, and D; the leakage-corrected ABCD equation gives S_A, and A minus S_A gives the fake-photon count.",
    "Second, ATLAS repeats the recoil analysis with Region-C candidates. Region C is isolated but non-tight and is dominated by neutral mesons and similar jet fragments, so its xJgamma distribution supplies the fake-tag recoil shape. Dividing by N_C makes that a per-candidate shape, and multiplying by N_bkg^A fixes how much of it belongs in Region A.",
    "Subtracting the scaled Region-C template from H_A removes the jets correlated with fake photon tags. In heavy-ion data, combinatoric jets are handled as a separate background before unfolding.",
    "[Sources]",
    "ATLAS Collaboration, Phys. Lett. B 789 (2019) 167, Sections 5.1 and 5.2.",
    "Original slide 7, PPG19_weeklyUpdate_8_14_26, revision clXJl_t1PNydnQ.",
  ]);
  return { slide, nodes };
}

async function main() {
  const outputDir = path.resolve(process.argv[2] || ".");
  await fs.mkdir(outputDir, { recursive: true });
  const stem = "slide07_purity_correction_method_flowchart_candidate";

  const math = await renderMath();
  const presentation = Presentation.create({ slideSize: { width: SLIDE_W, height: SLIDE_H } });
  const { slide, nodes } = buildSlide(presentation, math);

  const pngPath = path.join(outputDir, `${stem}.png`);
  const pptxPath = path.join(outputDir, `${stem}.pptx`);
  const layoutPath = path.join(outputDir, `${stem}_layout_nodes.json`);
  const artifactLayoutPath = path.join(outputDir, `${stem}_artifact_layout.json`);
  const manifestPath = path.join(outputDir, `${stem}_manifest.json`);
  const speakerPath = path.join(outputDir, `${stem}_speaker_script.md`);

  const preview = await presentation.export({ slide, format: "png", scale: SCALE });
  await saveBlob(preview, pngPath);
  const pythonBin = process.env.PYTHON_BIN || "python3";
  execFileSync(pythonBin, [
    "-c", "from PIL import Image; import sys; p=sys.argv[1]; Image.open(p).convert('RGB').save(p)", pngPath,
  ], { stdio: "pipe" });

  const artifactLayout = await presentation.export({ slide, format: "layout" });
  await saveBlob(artifactLayout, artifactLayoutPath);
  const pptx = await PresentationFile.exportPptx(presentation);
  await pptx.save(pptxPath);

  const layoutPayload = {
    slide_size: [SLIDE_W * SCALE, SLIDE_H * SCALE],
    title_axis_x: 68,
    minimum_audience_font_px: 38,
    minimum_title_font_px: 80,
    nodes,
  };
  await fs.writeFile(layoutPath, `${JSON.stringify(layoutPayload, null, 2)}\n`, "utf8");

  const speaker = [
    "# Purity Correction Method Choice — speaker script",
    "",
    "ATLAS treats the purity correction as two linked subtractions.",
    "",
    "First, I count candidates in the standard tight or non-tight and isolated or non-isolated ABCD regions. Region A is the tight isolated signal region. Signal simulation tells me how often genuine photons leak from A into B, C, and D. Solving the leakage-corrected ABCD equation gives S_A, the genuine prompt-photon count in A, and A minus S_A gives the fake-photon count.",
    "",
    "Second, I remove the recoil jets associated with those fake tags. H_A contains both genuine photon-jet events and fake-photon dijet events. Region C is isolated but non-tight and is dominated by neutral mesons and similar jet fragments, so H_C divided by N_C supplies the per-candidate fake-tag recoil shape. I multiply that shape by the fake-photon count in A and subtract it bin by bin from H_A.",
    "",
    "The result is the purity-corrected reconstructed xJgamma distribution. In heavy-ion data, the combinatoric-jet subtraction remains a separate correction before unfolding.",
    "",
  ].join("\n");
  await fs.writeFile(speakerPath, speaker, "utf8");

  const manifest = {
    artifact_id: "the219_slide07_purity_correction_method_flowchart_candidate",
    title: "Purity Correction Method Choice",
    source_deck: "PPG19_weeklyUpdate_8_14_26",
    source_presentation_id: "19pjXOc1CJ1Ed7xtw2ZrXffztpWkzCz4lN13L7CHPdQE",
    source_slide_object_id: "g3f6f210f969_0_365",
    source_revision_id: "clXJl_t1PNydnQ",
    output_state: "local candidate only; Google Slides not mutated",
    content_contract: [
      "All visible method content from source slide 7 is preserved.",
      "Title remains exactly 'Purity Correction Method Choice'.",
      "Body text is 14-16 pt equivalent on a 1280x720 slide canvas.",
      "Stage 1 estimates S_A and N_bkg^A with leakage-corrected ABCD.",
      "Stage 2 uses Region C for the fake-tag xJgamma shape, scales it by N_bkg^A, and subtracts it from Region A.",
    ],
    source_claims: [
      "ATLAS Collaboration, Phys. Lett. B 789 (2019) 167, Sections 5.1-5.2.",
      "PPG12 leakage-corrected ABCD method as represented in the current analysis infrastructure.",
    ],
    outputs: { pngPath, pptxPath, layoutPath, artifactLayoutPath, speakerPath },
  };
  await fs.writeFile(manifestPath, `${JSON.stringify(manifest, null, 2)}\n`, "utf8");

  console.log(JSON.stringify({ pngPath, pptxPath, layoutPath, artifactLayoutPath, manifestPath, speakerPath }, null, 2));
}

main().catch((error) => {
  console.error(error.stack || String(error));
  process.exit(1);
});
