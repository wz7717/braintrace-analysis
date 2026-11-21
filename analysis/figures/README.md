# Standalone manuscript figures

The seven programs use the saved aggregate tables or, for Figure 1, the reference schematic. Matplotlib's bundled DejaVu fonts render the Figure 4 title and the Figure 7A input label in PNG output. Both scripts retain native text in vector output. The [v38 render receipt](FIGURES_4_7_V38.json) identifies the scripts, inputs, fonts and PNGs checked against both v38 Word documents.

The [Figure 1 receipt](FIGURE_1_ACCEPTANCE_20261120.json) identifies the Ensembl 97 schematic. The [17 November figure record](FIGURE_ACCEPTANCE.json) and [20 November Figure 4 and 7 record](FIGURES_4_7_ACCEPTANCE_20261120.json) describe their respective editions. [Code origins](PLOT_CODE_PROVENANCE.csv) · [Input provenance](PLOT_INPUT_PROVENANCE.csv)

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
