# case_studies/keloid

Keloid fibroblast application of `represolve` (manuscript Results,
"Calibrated keloid prioritization"). The analysis scripts live at the
repository root of `scar_benchmark/` (they orchestrate large public
datasets, not part of the installable package):

- `j17_multicohort_signature_screen.py` — 21-donor, 4-cohort keloid
  signature + Hs27/RPE1 CRISPRa screen.
- `j25_keloid_null_calibration.py` / `j32_keloid_clean_signature.py` —
  fixed-size permutation nulls and signature-variant arms.
- `j37_hs27_guidesplit_repeats.csv` — guide-split replicate agreement
  (near-floor blind check).
- `65_r2_donor_sensitivity.py` — donor-composition sensitivity
  (scar-control exclusion, minimum-cell threshold).
- `j91_headline_tf_repcen.csv` — per-candidate guide-split rep_cen.

The package itself is application-agnostic; nothing under `src/` depends
on keloid data.
