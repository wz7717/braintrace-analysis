# BrainTrace analysis and reproducibility materials

The November 20 GitHub edition includes the [Ensembl 97 source refit, external-cohort results and acquisition notes](analysis/research/macaque_ensembl97/README.md). The [web application at the matching commit](https://github.com/wz7717/braintrace-web/tree/45666e9dca04df29a1541c88ba8aaff36a0213dc) uses the same model, lower references and input mappings. The [November 17 archive](https://doi.org/10.5281/zenodo.24311372) preserves the preceding analysis package.

Figures 1, 4 and 7 have checked-in drawing scripts for the current manuscript. Figure 7 uses the revised source and external aggregate inputs. [Plotting instructions and render receipts](analysis/figures/README.md)

The [paper version binding](PAPER_VERSION_BINDING.json) identifies the manuscript, supplements and application source for this Ensembl 97 edition. See the [paper and evidence index](PAPER_REPRODUCTION_MAP.csv) for its supporting analyses.

| Material | Use |
| --- | --- |
| [Software](software/README.md) | Frozen macaque Network/Group/Exact inference with 120 references and two prepared GEO examples. The existing human CEL release has its own frozen asset list. |
| [Study analysis](analysis/README.md) | S31 and S32 research entries, the B17 real aggregate and synthetic example, S33 code and evaluation summaries, environment locks and authorized-source instructions. |
| [AHBA target inputs](analysis/ahba/README.md) | Byte-identical post-transform target array and a minimal index. The fixed public model supplies feature and class order. The donor-equal Bo reference centroids required by Pearson have a separate access condition. |
| [Figures 1–7](analysis/figures/README.md) | Seven independent plotting scripts and nine saved aggregate input tables. The historical seven-image check used the 15 November 2026 manuscript. The November 20 receipts bind Figures 1, 4 and 7. |
| [Data access and licences](DATA_ACCESS_AND_LICENCES.md) | Fixed source links and SHA-256, object-level terms and restrictions on exact study inputs. |

The starting point for a figure redraw is its saved numerical table or schematic. The availability and exact input requirements for a study replay are in the route notes. The authors' historical replay receipts record their specific starting point and environment.

[Integration record](INTEGRATION_PROVENANCE.json) · [Public file manifest](PUBLIC_FILE_MANIFEST.csv) · [Content SHA-256](SHA256SUMS.txt) · [Released software demo](https://braintrace-web.streamlit.app/)
