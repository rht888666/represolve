# Figure/result → code mapping

Every quantitative claim in the manuscript maps to a script and an output
file in the parent repository (`scar_benchmark/`). Run order and paths are
frozen in `IMETA_REBUILD_20261118.md` (evidence ledger).

## Benchmark half

| Result | Script | Output |
|---|---|---|
| Hs27 cross-guide benchmark (rep_cen, lw/svd/moppca) | `j11_crossguide_benchmark.py` | `results/j11_crossguide_benchmark.csv`, `results/j11.log` |
| RPE1 cross-guide benchmark | `j11` RPE1 variant | `results/j11_crossguide_benchmark_rpe1.csv` |
| K562 essential fold-frozen + nn baselines | `j22_k562_foldfrozen_centroid.py` | `results/j22_k562_foldfrozen_centroid.csv` |
| K562 GWPS genome-scale benchmark (9,211 perturbations) | `j22` on `replogle_k562_gwps_halves.npz` | `results/j22_k562_gwps_foldfrozen_centroid.csv` |
| Gem-partition sensitivity | `j24_build_split_halves.py`, `j23`/`j26` | `results/j23_k562_split_norm_summary.csv`, `results/j26_k562_logcpm_summary.csv` |
| Mask sensitivity (4 masks) | `j31` | `results/j31_k562_mask_sensitivity.csv` |
| Perturbation-bootstrap CIs | `j30_k562_bootstrap_ci.py` | `results/j30_k562_bootstrap_ci.csv` |
| Decile ranking (essential + GWPS + logCPM) | `j33`, `j38` | `results/j33_k562_decile_ranking*.csv`, `results/j38_gwps_decile_ranking.csv` |
| In-sample k selection | `j36_k562_insample_k.py` | `results/j36_k562_insample_k*.csv` |
| 10× guide-split repeats (Hs27) | `j37_hs27_guidesplit_repeats.py` | `results/j37_hs27_guidesplit_repeats.csv` |
| GEARS honest-baseline audit | `j34_gears_trainmean_baseline.py` | `results/gears/j34_gears_trainmean_eval.csv` |
| Nonlinear (DAE) completion baseline, essential + GWPS | `j47_k562_dae_foldfrozen.py` | `results/j47_k562_dae_foldfrozen.csv`, `results/j47_gwps_dae_foldfrozen.csv` |
| External halves builder (scPerturb datasets) | `j48_build_external_halves.py` + VM `build_norman.py` | `norman_singles_gemhalves.npz` |
| Norman external replication benchmark | `j49_external_benchmark.py` | `results/j49_norman_foldfrozen.csv` |
| RPE1 (Replogle) external benchmark | `j49_external_benchmark.py` on `rpe1_batchhalves.npz` | `results/j49_rpe1_foldfrozen.csv` |
| sci-Plex3 chemical-modality benchmark | `j49_external_benchmark.py` on `sciplex3_rep_halves.npz` (rep1/rep2 split) | `results/j49_sciplex3.csv` |
| Cross-dataset model-ranking table | merge of j22/j47/j49 | `results/j51_crossdataset_model_ranking.csv` |
| Pre-declarable regime stratification (obs_agreement) | `j50_observed_gene_agreement.py` | `results/j50_obs_agreement.csv` |
| CPA unseen-combination audit (Norman doubles) | `cloud/cpa_norman.py` on `kvc-gpu2` | `results/checkpoints/j05_cpa_norman/cpa_norman_metrics.csv` (GCS) |
| DAE on Norman/RPE1 | `j47` variants | `results/j47_norman_dae.csv`, `results/j47_rpe1_dae.csv` |

## Keloid application half

| Result | Script | Output |
|---|---|---|
| 21-donor 4-cohort signature | `13_export_stage5b_bundle.py` + pipeline | `cloud_stage5b/keloid_21donor_*` |
| CRISPRa signature screen + permutation null | `exp4/exp5` scripts | `results/exp4_crispra_keloid_screen.csv`, `results/exp5_*` |
| DoRothEA regulon enrichment | `j40` | `results/j40_dorothea_regulon_enrichment.csv` |
| Open Targets associations | `j41` | `results/j41_opentargets_keloid_associations.csv` |
| L1000CDS2 reversal + random-signature null | `j42`, `j44` | `results/j42_l1000cds2_*`, `results/j44_l1000_randdir_null.csv` |
| Extended all-pairs coloc + fine-mapping | `43_full_allpairs_coloc_and_finemapping.py` | `results/j43_coloc_extension_candidates.csv` |
| Evidence matrix | `j45` | `results/j45_candidate_evidence_matrix.csv` |
| Candidate fibroblast expression (HPA substitute) | `j52_candidate_fibroblast_expression.py` | `results/j52_candidate_fibroblast_expression.csv` |
| Type-I calibration of screen null | `j46` | `results/j46_null_type1_check.csv` |

## AlphaGenome 2 half

| Result | Script | Output |
|---|---|---|
| Credible-set construction | `ag2_02_build_credible_sets.py` | `results/alphagenome2/ag2_02_credible_sets.csv` |
| Variant→gene nominations | `ag2_11` | `results/ag2_11_v2g_nominations.csv` |
| GTEx eQTL sign calibration | `ag2_12` | `results/ag2_12_gtex_sign_calibration.csv` |
| Donor-level conditioned stats | `ag2_13` | `results/ag2_13_donor_level_stats.csv` |
| Matched empirical null | `ag2_14` + `ag2_14_null_variants.py` | `results/ag2_14_empirical_null_{RNA_SEQ,DNASE}.csv` |
| Perturbation-state conditioning + bridge | `ag2_09`, `ag2_10` | `results/alphagenome2/ag2_09_*.csv`, `results/ag2_10_*` |

## Manuscripts and provenance

- Working outline/evidence ledger: `IMETA_REBUILD_20261118.md`
- Analysis audit trail: `AUDIT_LOG.md`, `BENCHMARK_LOG.md`
- Cloud execution + checkpoint policy: `AGENTS.md`
