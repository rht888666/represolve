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
| Fibrosis/ECM-module polarity (6 TF hits) | `j99_tf_fibrosis_module.py` | `results/j99_tf_fibrosis_summary.csv` |
| Package v1.0.0 release prep | `represolve/pyproject.toml`, `README.md`, `CHANGELOG.md` | `submission/RELEASE_AND_PUBLISH.md` |
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

## Additions 2026-11-18 (round-1 review fixes)

| Result | Script | Output |
|---|---|---|
| Pairwise ranking bootstrap CIs | `j53_pairwise_ranking_ci.py` | `results/j53_pairwise_ranking_ci.csv` |
| obs_agreement proxy, all contexts | `j50_observed_gene_agreement.py` | `results/j50_obs_agreement{,_norman,_rpe1,_gwps,_sciplex3}.csv` |
| CPA OOD-combo reconstruction check | `cloud/cpa_norman_ood_recon.py` | `results/cpa_norman_ood_recon.csv` |
| Frangieh guide-split context | `build_frangieh_guidehalves.py` on scPerturb `FrangiehIzar2021_RNA.h5ad` | `results/j73_frangieh_norm_benchmark.csv`, `j75_dae_frangieh_guidehalves.csv`, `j78_frangieh_merged.csv` |
| Tian CRISPRi batch-split context | inline build on scPerturb `TianKampmann2021_CRISPRi.h5ad` (batch odd/even, log1p CPM10k) | `results/j74_tian2021i_norm_benchmark.csv`, `j75_dae_tian2021i_oddeven_halves.csv`, `j78_tian2021i_merged.csv` |
| RPE1 E150 gem-split context | inline build on `RPE1_E150_all_genes.h5ad` | `results/j72_rpe1e150_benchmark.csv`, `j75_dae_rpe1_e150_gemhalves.csv`, `j78_rpe1e150_merged.csv` |
| Tahoe-100M plate3 sample-split context | inline build on `tahoe100m/plate3.h5ad` | `results/j77_tahoe_benchmark.csv`, `j75_dae_tahoe.csv`, `j78_tahoe_merged.csv` |
| scGPT zero-shot audit | `j76_scgpt_zeroshot.py` (fold 0, binned inputs) | `results/j76_scgpt_norman_f0.csv`, `results/j76_scgpt_k562_f0.csv` |
| obs_agreement proxies, new contexts | `j50_observed_gene_agreement.py` | `results/j50_obs_agr_*.csv` |

## Additions 2026-11-20 (round-2 review fixes)

| Result | Script | Output |
|---|---|---|
| Cross-half training ablation (merged vs half-A fit) | `62_ablations_resume.py` (lean `61_ablations_lean.py`) | `results/j91_crosshalf.csv` |
| Mask=0.2 sensitivity, all contexts | `62_ablations_resume.py` | `results/j91_mask02.csv` |
| sci-Plex3 drug-GroupKFold | `62_ablations_resume.py` | `results/j91_sciplex3_groupkfold.csv` |
| Perturbation-density metric (nn10 distance) | `62_ablations_resume.py` | `results/j91_density.csv` |
| sci-Plex3 N=145 downsample (masked-MLP audit) | `61_ablations_lean.py` | `results/j91_sciplex3_mlp_N145.csv` |
| Keloid donor metadata + Mann-Whitney | inline (`j28_signature_external_scores.csv`) | `results/j91_donor_metadata_table.csv` |
| Leave-one-cohort-out CV AUC | `62` LOCO block | `results/j91_loco_cv_auc.csv` |
| Headline-TF guide-split rep_cen (10 repeats) | inline on `hs27` gem halves | `results/j91_headline_tf_repcen.csv` |
| obs_agreement end-to-end routing regret | inline (Table 3 rules vs oracle) | `results/j91_routing_regret.csv` |
| scGPT zero-shot 5-fold audit, 3 contexts | `cloud/j92_scgpt_5fold.py` on `kvc-gpu2` | `results/checkpoints/j92_scgpt_5fold/j92_summary.csv` + per-context CSVs |

## Additions 2026-11-20 (external validation, fine-tuned scGPT, visuals)

| Result | Script | Output |
|---|---|---|
| Out-of-sample halves (HepG2, Jurkat, sci-Plex4, Tian CRISPRa) | `70_build_external_halves.py` | `zenodo_upload/data_external/*_halves.npz` |
| Pre-registered routing predictions (SHA-256 + UTC) | `71_preregister_routing.py` | `results/j94_preregistered_routing.csv` |
| Out-of-sample hidden-gene scoring & regret | `72_external_benchmark.py` | `results/j95_external_validation.csv`, `results/j95_*_perpert.csv` |
| End-to-end fine-tuned scGPT 5-fold audit (3 contexts) | `cloud/j96_scgpt_finetune.py` on `kvc-j96` (L4) | `results/checkpoints/j96_scgpt_finetune/j96_summary.csv` + per-context CSVs |
| Redesigned Figure 1A/1B & 2D phase-map Graphical Abstract | `represolve/notebooks/fig1_protocol.ipynb`, `make_graphical_abstract.py` | `figures/fig1a_schematic.png`, `figures/fig1b_same_vs_indep.png`, `figures/graphical_abstract.png` |
| Publication-grade Figures 2–6 redesign (incl. Fig 3C scGPT fine-tuned & Fig 6A/B 9+4 phase map + external validation) | `make_figures_v2.py` | `figures/fig2..6*.png`, `figures/_comp_fig2..6.png` |
| Unseen-perturbation replicate-resolved generalization audit (4 genetic screens) | `73_unseen_perturbation_audit.py` | `results/j98_unseen_perturbation_audit.csv`, `results/j98_unseen_*_perpert.csv` |
| Interactive video demo + 5 presentation slides + bilingual SRT | `represolve/demo/run_video_demo.py`, `make_video_slides.py` | `zenodo_upload/VIDEO_SCRIPT.md`, `zenodo_upload/video_slides/*` |
