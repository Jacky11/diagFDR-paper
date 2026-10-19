# Figures

Figures 1, 3, 4 and 5 of the paper are composites assembled from the panels
written by the scripts of `real_data/` (one PDF per panel). Two steps are scripted
here.

## Figure 2 (workflow and report)

1. `Rscript figure2B_report.R` (needs `diagFDR` >= 0.1.2 and
   `../real_data/out/universes.rds`) runs `dfdr_run_all()` on three lists of the
   basic DIA-NN search (global, min_r, run x precursor), marked as two-pass
   exports, prints the headline values (D = 832 / 2001 / 3751; IPE 0.075 /
   0.075 / 6e-5) and renders the HTML report to `out/report/`.
2. `python figure2B_screenshots.py` (Playwright + Chromium) captures the title,
   table of contents, executive summary and headline table of that report.
      Crop the executive summary above "Quick Metrics Dashboard" (1120 px at
   scale 2), the headline table above "Click to show" (600 px) and the table of
   contents to its text width.
3. `python html_to_pdf.py $(pwd)/figure2.html Figure2.pdf 1600 1600` lays out
      panel A (workflow) and panel B (report captures) and prints the PDF.
4. `node figure2_pptx.js` (needs `pptxgenjs`, and the four captures in the
   same folder) writes `Figure2_editable.pptx`, an editable PowerPoint version:
   panel A is built from native shapes and text boxes, panel B places the
   report captures as pictures with native arrows and labels.

## Figures 3 to 5: labels

`relabel_figures_3_to_5.py` replaces, in the PowerPoint-exported composites,
the labels that called the decoy-support index a "CV" (Figure 3D legend and
subtitle, Figure 4A title, subtitle and right axis, Figure 5B title and axis)
and the y-axis of Figure 4D (`D_b / T_b`, the comparator of the IPE). It redacts the
original text under each box (PyMuPDF), so the old labels are removed from the
text layer and not only hidden, then draws the new text from a reportlab overlay
(`relabel_overlay.py`, needs `pymupdf`, `reportlab` and the Liberation Sans and
DejaVu Sans fonts). The panels themselves are unchanged.

Inputs expected by the scripts: `F3_orig.pdf`, `F4_orig.pdf`, `F5_orig.pdf`
(the composites before relabelling: scope audit, calibration and stability,
MS2Rescore), and for `figure2.html` the captures `title.png`, `toc_crop.png`,
`exec_crop.png`, `headline_crop.png` in the same folder.
