# BrainTrace analysis and reproducibility materials

The [Zenodo archive (version DOI 10.5281/zenodo.24311372)](https://doi.org/10.5281/zenodo.24311372) provides the manuscript's analysis package, `BrainTrace_Analysis_c5ac7bc_Prepublication_20261117.zip`, and the accompanying version record.

To redraw Figure 4 in the current manuscript, use `analysis/figures/plot_figure_4.py` from the extracted Zenodo ZIP. The [script at GitHub commit `c5ac7bc`](https://github.com/wz7717/braintrace-analysis/blob/c5ac7bcff8d42c91942ef1f1207daf8f74642eb0/analysis/figures/plot_figure_4.py) produces the earlier version of Figure 4.

The [paper version binding](PAPER_VERSION_BINDING.json) identifies the manuscript, supplement and deployed software for the 15 November 2026 version. See the [paper and evidence index](PAPER_REPRODUCTION_MAP.csv) for its supporting analyses.

| Material | Use |
| --- | --- |
| [Software](software/README_CANDIDATE.md) | Frozen macaque Network/Group/Exact inference with 120 references and two prepared GEO examples. The existing human CEL release has its own frozen asset list. |
| [Study analysis](analysis/README.md) | S31 and S32 research entries, the B17 real aggregate and synthetic example, S33 code and evaluation summaries, environment locks and authorized-source instructions. |
| [AHBA target inputs](analysis/ahba/README.md) | Byte-identical post-transform target array and a minimal index. The fixed public model supplies feature and class order. The donor-equal Bo reference centroids required by Pearson have a separate access condition. |
| [Figures 1–7](analysis/figures/README.md) | Seven independent plotting scripts and nine saved aggregate input tables. Their PNG outputs matched the seven figures in the 15 November 2026 manuscript under the recorded environments. |
| [Data access and licences](DATA_ACCESS_AND_LICENCES.md) | Fixed source links and SHA-256, object-level terms and restrictions on exact study inputs. |

The starting point for a figure redraw is its saved numerical table or schematic. The availability and exact input requirements for a study replay are in the route notes. The authors' historical replay receipts record their specific starting point and environment.

[Integration record](INTEGRATION_PROVENANCE.json) · [Public file manifest](PUBLIC_FILE_MANIFEST.csv) · [Content SHA-256](SHA256SUMS.txt) · [Released software demo](https://braintrace-web.streamlit.app/)
