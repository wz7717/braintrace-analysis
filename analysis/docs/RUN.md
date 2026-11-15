# Running the included analyses and figures

Paths below are relative to the unpacked analysis directory. Use a separate output directory and verify any acquired frozen files against the recorded locks. [Data acquisition and access](DATA_ACCESS.md) · [Fixed analysis settings](../environment/ANALYSIS_SETTINGS.json)

## B17 synthetic Pearson check

On Windows x64 with CPython 3.11.9, create an environment from the supplied lock and run the synthetic example.

~~~powershell
py -3.11 -m venv .venv
& ./.venv/Scripts/python.exe -m pip install --require-hashes -r environment/b17-synthetic-windows-py311.lock
& ./.venv/Scripts/python.exe -B examples/b17_pearson_smoke.py
& ./.venv/Scripts/python.exe -B code/b17/run_comparison.py --help
~~~

The expected nine scores are 1 followed by eight values of -0.125 within absolute tolerance 1e-12. The historically frozen adult GSE25219 rows appear in the B17 results CSV. Real scoring requires the 46 frozen PLM queries and final source references. The B17 path-only adapter selects the original R 4.6.1 and SingleR 2.14.2 environment using BRAINTRACE_B17_RSCRIPT and BRAINTRACE_B17_R_LIBRARY. The R 4.6.2 lock is for S31. [B17 six rows](../results/B17_SUMMARY.csv)

## Study routes

| Route | Entry with authorized inputs | Evidence to check |
| --- | --- | --- |
| S31 SingleR | Run research/human/run_comparison.py with --config and --output. It uses the paired research/human/score_singler.R from the same code directory. An optional --work-dir may use a short ASCII path. | Locked seven objects, 1,856 queries, 200 ordered markers and 16,704 scores. Original quantile 0.8, no fine-tuning or pruning. |
| S31 UKBEC | Run research/human/ukbec_rank_rf500_diagnostic.py with --predictions and --output against the original saved five-fold rows. | Check all 1,083 saved predictions and the original RF500 diagnostic strata. |
| S32 Bo source | From research, use python -B -m internal.run_fixed_gate --config CONFIG.json --out OUTPUT_DIR. Use a new output directory and the original frozen source/fold inputs. | Nine frozen outer folds, 819 records and the original in-pool, out-of-pool and unsupported strata. |
| S32 GEO / PsychENCODE | From research, use python -B -m external.shared_head_external or python -B -m external.shared_head_psychencode, with --config and --output. Reuse the matching head, masks and labels. | Two GEO cohorts with 675 comparable, 268 unscored. Adult PsychENCODE 34 comparable and 43 unscored at the fine tier. |
| S32 AHBA | From research, external.matched_ahba_network uses the official inputs, prepared metadata, source B5 members, model and exact run config. The separate split target stores the already-transformed matrix and index. | The historical post-transform run reproduced 2,310 logits and the three result CSVs. This candidate omits Bo centroids and the row-reporting full library total from the split index. [AHBA split](../ahba/README.md) |
| S33 ACT | Run act/replay_s33.py with --config and --output. Use the held minimal roster, exact TPM bodies, fixed model and the matching rank background. The exact historical comparator also needs the authorized original directory. | Nine classes for 21 specimens, 189 logits and 189 scores with original rank and status. Table S50 reports 11 donor-equal. [ACT access and execution](../act/README.md) |

## Environment

- [S31/S32 Windows Python lock](../environment/s31-s32-windows-py311.lock) uses Windows x64 CPython 3.11.9. Use the thread controls specified by the original run.
- [S31 compatible R download lock](../environment/S31_R_PUBLIC_COMPATIBLE.csv) identifies R 4.6.2 UCRT, Bioconductor 3.23 and the exact compatible archives. Verify SHA-256 before installing into an isolated R_LIBS_USER with repos = NULL. The corresponding three public Rcpp, generics and matrixStats binaries differ at the byte level from the older same-version builds.
- [S33 Windows Python lock](../environment/s33-windows-py311.lock) identifies the eight wheels from the verified later replay. The original run did not save an independent helper SHA. [Original and compatible environments](../environment/S31_S33_ENVIRONMENTS.json)

## Figures 1–7

The seven independent programs render from the saved source values. Use the separate instructions and exact rendering environment of the relevant figure. [Figure runbook](../figures/README.md)
