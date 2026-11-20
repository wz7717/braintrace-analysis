# Macaque Ensembl 97 source refit and fixed transfer

The November 20 revision evaluates 819 Bo et al. specimens across nine held-out donors. The source annotation reconstructs the M. fascicularis 5.0 Ensembl release 97 gene map. All 28,415 stable IDs in the count matrix match. Within every training split the detection threshold is at least 10% of training specimens. Network feature selection and model fitting use the eligible training genes. The deployable model and the 120 unordered Network Top3 references use all nine source donors.

## Source run

Use Python 3.11 with NumPy 1.26.4, pandas 2.2.3, SciPy 1.13.1, scikit-learn 1.5.2, joblib, threadpoolctl and openpyxl. From this directory, run the following command. Supply the published Bo count matrix and the specimen workbook under their original filenames in `external_data/Bo2023/` below the input directory. The `source_data.py` reader specifies their filenames and workbook sheet.

```bash
python source_analysis/run_source.py --science-root /path/to/input --output source_analysis/results --workers 2
```

The program fits its selection grid and produces the held-out predictions and final source assets in a new output directory. The final scaling parameters and linear classifier, gene axes, lower-tier references, platform mappings and two runnable GEO examples are also distributed under [software/runtime_assets](../../../software/runtime_assets). See [annotation provenance and upstream sources](source_analysis/MAPPING_SOURCE_AND_CHECKS.json), the [Bo project](https://github.com/WangLab-SINH/Macaque_Brain_Transcriptome_MRI) and the [aggregate source results](aggregates/SOURCE_REFIT.json).

## Fixed inference and study replay

For an individual compatible input or the two GEO examples, use the [application](../../../software/README.md). For a study-level replay, set `BRAINTRACE_STUDY_ROOT` to the directory that holds the historical `BrainTrace_algorithm_research_20261015` and associated cohort inputs. The aggregate files below can be read directly.

```bash
python external_analysis/transfer.py
python external_analysis/endpoints.py
python external_analysis/retrieval.py
python scientific_results/source_tables.py
python scientific_results/source_donor_details.py
```

Study replay uses the original GEO and adult PsychENCODE count matrices, platform gene correspondence, historical truth and evaluation rosters, resampling draws and baseline predictions. The two-donor AHBA route also uses the original Allen archives and the frozen coverage roster with paths resolved for the local installation. Source references come from the preceding source run. The optional Group stability and annotation-comparison scripts use prior study outputs as documented in their path declarations. The Bo authors' original files, recipient-specific rosters and the historical expression matrices are obtained separately.

## Numerical outputs

* [Source three-tier rates](aggregates/SOURCE_TIER_RATES.csv) and [external endpoints](aggregates/FROZEN_ENDPOINTS.csv) distinguish the original and revised routes. In the external table, `corrected_*` denotes the Ensembl 97 result.
* [Source equal-budget retrieval](aggregates/TABLE_7_S38.csv) and [external equal-budget retrieval](aggregates/RETRIEVAL_MATCHED.csv) contain the numerical results for Table 7 and S38. The ambiguous GEO rows use each cohort’s own donor and candidate counts.
* [Observed mapping coverage](aggregates/coverage.csv) reports the 5,000-feature Network axis and 21,668-feature lower axis. Each GEO matrix and the PsychENCODE common annotation have their own count background.
* [Figure 7](../../figures/README.md) can be drawn directly from the two distributed aggregate tables. [Figure S10 points](aggregates/FIGURE_S10_POINTS.csv) and the other [aggregate tables](aggregates/) are also supplied.

The historical releases and their notices remain at their original Git revisions. Third-party data, annotations and software retain their source terms.
