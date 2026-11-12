# Provenance of the results in this folder (11 November 2026, after the complete R run)

| Item | Source |
|---|---|
| tables/TableS1-S14c (all `Table*.csv` except S15-S17), tables/manuscript_numbers.csv, figures/*, data/processed/*.rds, logs/sessionInfo.txt | complete run of `run_all.R` in R 4.6.0 on the author's machine (11 November 2026: steps 01-06 at 13:34-13:57, steps 07-10 at 15:12-15:19, step 11 at 15:19-15:21) |
| tables/TableS16 (GIGASTROKE rows), TableS17 (cardioembolic rows) | `11_reverse_and_af_mvmr_token.qmd` (R); written into the tables by `tools/recompute_s17.py` and from `data/processed/reverse_mr.rds` |
| tables/TableS15*, TableS16 (MEGASTROKE rows), TableS17 (other outcomes), TableS17b, TableS17d, manuscript_numbers_extra.csv | Python tools (`tools/ld_proxy_sensitivity.py`, `tools/reverse_and_af_mvmr_api.py`, `tools/mr_api_lib.py`), validated against R (Table S17c) |
| tables_baseline_AF130/ | AF tables and number list of the v9 run (130-lead cache), kept as evidence |
| tables/independent_*_results.csv, af_update_validation_vs_v9.csv | Python checks of 10 November |
