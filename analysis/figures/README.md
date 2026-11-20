# Standalone manuscript figures

For Figure 1, `plot_figure_1.py` draws the 20 November 2026 schematic with GTEx, HPA and Allen ACT. The [Figure 1 render receipt](FIGURE_1_ACCEPTANCE_20261120.json) records the exact PNG. The earlier seven-figure acceptance record refers to its original manuscript.

The original seven-image check used the stored figure values and schematic for the 15 November manuscript. The receipt above binds the revised Figure 1. The [November 20 Figure 4 and Figure 7 receipt](FIGURES_4_7_ACCEPTANCE_20261120.json) binds the public drawing scripts to the corresponding embedded images. Figure 4 uses the archived November 17 script. Figure 7 uses both revised aggregate inputs. Each program runs separately. The input tables hold the published aggregate task strata, counts and saved intervals relevant to that figure. [Code origins](PLOT_CODE_PROVENANCE.csv) · [Input provenance](PLOT_INPUT_PROVENANCE.csv) · [Seven-image QA receipt](FIGURE_ACCEPTANCE.json)

| Figures | Checked rendering environment |
| --- | --- |
| 1 | Windows x64 CPython 3.12.14, ReportLab 4.4.9, Poppler pdftoppm 26.07.0, and local Arial regular and bold. The receipt gives the original font digests. |
| 2–4, 7 | Windows x64 CPython 3.11.9, Matplotlib 3.10.7. [Checked package versions](requirements-figures-2-4.txt) |
| 5–6 | Windows x64 CPython 3.12.14, Matplotlib 3.11.2. [Checked package versions](requirements-figures-5-7.txt) |

Run from this directory with the interpreter for the figure. On Windows, examples are

~~~powershell
# Use the Figure 1 environment with ReportLab and Poppler.
python -B plot_figure_1.py --font-regular "$env:WINDIR/Fonts/arial.ttf" --font-bold "$env:WINDIR/Fonts/arialbd.ttf" --pdftoppm pdftoppm --output plots/Figure_1.png

# Use the Figures 2–4 and 7 environment.
foreach ($figureNumber in @(2, 3, 4, 7)) {
    python -B "plot_figure_$figureNumber.py" --output "plots/Figure_$figureNumber.png"
}

# Use the Figures 5–6 environment.
foreach ($figureNumber in 5..6) {
    python -B "plot_figure_$figureNumber.py" --output "plots/Figure_$figureNumber.png"
}
~~~

Figures 2–7 use their bundled tables by default. Use --inputs to point to an equivalent copy. Figure 1 converts a temporary PDF to the PNG at 600 dpi. Different font or renderer builds can affect exact output bytes; the acceptance identity applies to the tested stack. The figure package contains script sources and tables, with PNG outputs generated into the chosen directory. The scholarly output identity is that of the embedded Word images under the manuscript hash in the receipt.
