// Editable PowerPoint version of Figure 2 (diagFDR workflow + report captures).
// Panel A is built from native shapes and text boxes; panel B uses the report
// captures as pictures, with native labels and arrows. Coordinates follow the
// 1600 x 1600 px layout of figure2.html (160 px = 1 inch).
const pptxgen = require("pptxgenjs");
const pres = new pptxgen();
pres.defineLayout({ name: "FIG2", width: 10, height: 10 });
pres.layout = "FIG2";
pres.title = "diagFDR Figure 2";
pres.author = "diagFDR authors";

const S = pres.ShapeType;
const px = (v) => v / 160;          // px -> inch
const pt = (v) => +(v * 0.45).toFixed(1); // px font size -> pt
const INK = "1B2333";
const SANS = "Arial";
const MONO = "Courier New";
const slide = pres.addSlide();
slide.background = { color: "FFFFFF" };

function box(name, x, y, w, h, fill, dashed) {
  slide.addShape(S.rect, {
    x: px(x), y: px(y), w: px(w), h: px(h),
    fill: { color: fill },
    line: { color: INK, width: 2.2, dashType: dashed ? "dash" : "solid" },
    objectName: name,
  });
}
function text(name, runs, x, y, w, h, opts = {}) {
  slide.addText(runs, Object.assign({
    x: px(x), y: px(y), w: px(w), h: px(h), margin: 0, valign: "top",
    fontFace: SANS, color: "111111", isTextBox: true, objectName: name,
  }, opts));
}
const r = (t, o = {}) => ({ text: t, options: o });
const br = (t, o = {}) => ({ text: t, options: Object.assign({ breakLine: true }, o) });
function badge(name, cx, cy, n) {
  slide.addShape(S.ellipse, {
    x: px(cx - 27), y: px(cy - 27), w: px(54), h: px(54),
    fill: { color: "1E9C43" }, line: { color: "1E9C43", width: 0 }, objectName: name,
  });
  text(name + " number", [r(String(n), { bold: true, color: "FFFFFF" })],
    cx - 27, cy - 27, 54, 54, { fontSize: pt(28), align: "center", valign: "middle" });
}
// cubic curve (from figure2.html) as an editable freeform line
function curve(name, p0, c1, c2, p1, arrow) {
  const x0 = Math.min(p0[0], c1[0], c2[0], p1[0]), y0 = Math.min(p0[1], c1[1], c2[1], p1[1]);
  const x1 = Math.max(p0[0], c1[0], c2[0], p1[0]), y1 = Math.max(p0[1], c1[1], c2[1], p1[1]);
  const rel = (p) => ({ x: px(p[0] - x0), y: px(p[1] - y0) });
  slide.addShape(S.custGeom, {
    x: px(x0), y: px(y0), w: Math.max(px(x1 - x0), 0.01), h: Math.max(px(y1 - y0), 0.01),
    line: { color: INK, width: 2.2, endArrowType: arrow ? "triangle" : undefined },
    points: [
      rel(p0),
      Object.assign(rel(p1), { curve: { type: "cubic", x1: px(c1[0] - x0), y1: px(c1[1] - y0), x2: px(c2[0] - x0), y2: px(c2[1] - y0) } }),
    ],
    objectName: name,
  });
}

// ---------------- Panel A ----------------
text("Label A", [r("A", { bold: true })], 16, 8, 60, 60, { fontSize: pt(46) });
text("Input title", [r("Input", { bold: true })], 120, 92, 210, 36, { fontSize: pt(28), align: "center" });
text("Output title", [r("Output", { bold: true })], 1215, 92, 375, 36, { fontSize: pt(28), align: "center" });

box("Input box", 10, 136, 440, 420, "EEF1FB");
text("Input text", [
  br("Any search engine output", { bold: true, fontSize: pt(22) }),
  r("Required: ", { bold: true }), br("id, is_decoy, score"),
  r("Used if present: ", { bold: true }), br("q-value, PEP, p-value, run"),
  r("External validation: entrapment or ground-truth labels", { color: "4A4A4A" }),
], 34, 152, 400, 160, { fontSize: pt(19), lineSpacingMultiple: 1.35 });
box("Adapters box", 34, 340, 392, 200, "EEF1FB", true);
text("Adapters text", [
  br("Provided adapters (Table S1)", { bold: true, fontFace: SANS, fontSize: pt(19), color: "333333" }),
  br('read_diann_parquet("report.parquet")'),
  br('read_spectronaut("report.tsv")'),
  br('read_dfdr_maxquant_msms("msms.txt")'),
  br('read_dfdr_mzid("file.mzid")'),
  r('read_mokapot_psms("t.psms.txt", "d.psms.txt")'),
], 46, 350, 372, 186, { fontFace: MONO, fontSize: pt(16.5), lineSpacingMultiple: 1.3 });

box("Install box", 520, 40, 640, 86, "E6E6E6");
text("Install text", [
  br("Installing and loading diagFDR in R", { bold: true, fontSize: pt(22) }),
  r('install.packages("diagFDR"); library(diagFDR)', { fontFace: MONO, fontSize: pt(16.5) }),
], 543, 54, 600, 64);

box("Step 1 box", 520, 140, 640, 92, "E8F8EE");
text("Step 1 text", [
  br("Import data (tool-agnostic)", { bold: true, fontSize: pt(22) }),
  r("Adapters or column mapping: ", { color: "4A4A4A", fontSize: pt(19) }),
  r("as_dfdr_tbl()", { fontFace: MONO, fontSize: pt(16.5), color: "4A4A4A" }),
], 543, 156, 560, 66);
badge("Step 1 badge", 1157, 186, 1);

box("Step 2 box", 520, 246, 640, 132, "E8F8EE");
text("Step 2 text", [
  r("Run diagnostics: ", { bold: true, fontSize: pt(22) }),
  br("dfdr_run_all()", { bold: true, fontFace: MONO, fontSize: pt(19) }),
  r("Scope (decoy ratio D/T per q-value column) • calibration (IPE, equal chance) • stability (D", { color: "4A4A4A" }),
  r("α", { color: "4A4A4A", subscript: true }),
  r(", FDP-CV benchmark, D", { color: "4A4A4A" }),
  r("α,win", { color: "4A4A4A", subscript: true }),
  r(", elasticity)", { color: "4A4A4A" }),
], 543, 262, 560, 106, { fontSize: pt(19), lineSpacingMultiple: 1.3 });
badge("Step 2 badge", 1157, 312, 2);

box("Optional box", 520, 392, 640, 76, "F6F7FB", true);
text("Optional text", [
  r("Optional: ", { bold: true }), r("pseudo-p-values from decoys  "),
  r("score_to_pvalue()", { fontFace: MONO, fontSize: pt(16.5) }),
], 543, 410, 600, 40, { fontSize: pt(19) });

box("Step 3 box", 520, 482, 640, 74, "E8F8EE");
text("Step 3 text", [r("Report and share", { bold: true })], 543, 498, 560, 34, { fontSize: pt(22) });
badge("Step 3 badge", 1157, 519, 3);

box("Bundle box", 1215, 136, 375, 200, "FDF3E4");
text("Bundle text", [
  br("Reproducible bundle", { bold: true, fontSize: pt(22) }),
  br("dfdr_write_report()", { fontFace: MONO, fontSize: pt(16.5) }),
  br("• tables (CSV), plots (PNG)"),
  r("• manifest + README + summary"),
], 1238, 152, 335, 170, { fontSize: pt(19), lineSpacingMultiple: 1.35 });

box("Report box", 1215, 356, 375, 200, "FDF3E4");
text("Report text", [
  br("HTML report", { bold: true, fontSize: pt(22) }),
  br("dfdr_render_report()", { fontFace: MONO, fontSize: pt(16.5) }),
  r("single entry point for readers and reviewers; alerts, not a verdict"),
], 1238, 372, 335, 170, { fontSize: pt(19), lineSpacingMultiple: 1.35 });

curve("Connector input-import", [450, 340], [490, 340], [480, 190], [520, 186], false);
curve("Connector report-bundle", [1184, 519], [1200, 519], [1195, 240], [1215, 236], false);
curve("Connector report-html", [1184, 519], [1200, 519], [1195, 460], [1215, 456], false);

// ---------------- Panel B ----------------
text("Label B", [r("B", { bold: true })], 16, 590, 60, 60, { fontSize: pt(46) });
const img = (name, file, x, y, w, wpx, hpx) =>
  slide.addImage({ path: file, x: px(x), y: px(y), w: px(w), h: px(w * hpx / wpx), altText: name, objectName: name });
img("Report title capture", "title.png", 20, 648, 430, 1820, 168);
img("Report table of contents capture", "toc_crop.png", 20, 706, 430, 994, 1098);
img("Report executive summary capture", "exec_crop.png", 620, 640, 960, 1820, 1120);
img("Report headline table capture", "headline_crop.png", 620, 1250, 960, 1820, 598);
slide.addShape(S.line, {
  x: px(455), y: px(760), w: px(155), h: 0,
  line: { color: INK, width: 2.2, endArrowType: "triangle" }, objectName: "Arrow to executive summary",
});
curve("Arrow to headline table", [455, 812], [560, 812], [520, 1325], [610, 1330], true);

pres.writeFile({ fileName: "Figure2_editable.pptx" }).then((f) => console.log("wrote", f));
