"""Frozen replicate-half datasets used in the manuscript.

Each entry describes one masked-completion benchmark context: two
independent measurement halves (S[0], S[1]) of the same perturbation
panel, stored as ``.npz`` with keys ``S`` (2, n_pert, n_genes),
``N`` (2, n_pert) cell counts, ``perts``, ``genes``.

Halves were built by ``j48_build_external_halves.py`` (scPerturb Zenodo
record 7041849) and the K562 builders described in
``paper_map/FIGURE_TO_CODE.md``. Checksums are in
``DATA_MANIFEST.csv`` next to this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class HalfDataset:
    name: str
    description: str
    split_axis: str          # what the two halves separate
    modality: str
    source: str


DATASETS = {
    "k562_essential": HalfDataset(
        name="k562_essential",
        description="Replogle Weismann 2022 K562 essential screen, "
                    "even/odd gem-group cell-pool halves",
        split_axis="gem_group (cell/lane pools; not donors)",
        modality="CRISPRi",
        source="Replogle et al. 2022, Cell; normalized single-cell h5ad",
    ),
    "k562_gwps": HalfDataset(
        name="k562_gwps",
        description="Replogle Weismann 2022 K562 genome-wide "
                    "perturb-seq, even/odd gem-group halves",
        split_axis="gem_group (cell/lane pools; not donors)",
        modality="CRISPRi genome-scale",
        source="Replogle et al. 2022, Cell",
    ),
    "norman": HalfDataset(
        name="norman",
        description="Norman et al. 2019 combinatorial screen, single-gene "
                    "perturbations, odd/even gem-group halves",
        split_axis="gemgroup",
        modality="CRISPRa",
        source="NormanWeissman2019_filtered.h5ad, scPerturb Zenodo 7041849",
    ),
    "rpe1": HalfDataset(
        name="rpe1",
        description="Replogle RPE1 screen, odd/even batch halves",
        split_axis="batch (56 batches, odd/even)",
        modality="CRISPRi",
        source="ReplogleWeissman2022_rpe1.h5ad, scPerturb Zenodo 7041849",
    ),
    "sciplex3": HalfDataset(
        name="sciplex3",
        description="Srivatsan et al. 2020 sci-Plex3 chemical screen, "
                    "rep1/rep2 biological-replicate halves across "
                    "drug x dose x cell line",
        split_axis="biological replicate (rep1 vs rep2)",
        modality="chemical perturbation",
        source="SrivatsanTrapnell2020_sciplex3.h5ad, scPerturb Zenodo 7041849",
    ),
    "frangieh": HalfDataset(
        name="frangieh",
        description="Frangieh et al. 2021 melanoma CRISPR screen, "
                    "guide-level halves",
        split_axis="guide split",
        modality="CRISPR",
        source="Frangieh et al. 2021, Nature Genetics; scPerturb Zenodo 7041849",
    ),
    "tian2021i": HalfDataset(
        name="tian2021i",
        description="Tian et al. 2021 CRISPRi screen, odd/even halves",
        split_axis="odd/even split",
        modality="CRISPRi",
        source="Tian et al. 2021; scPerturb Zenodo 7041849",
    ),
    "rpe1_e150": HalfDataset(
        name="rpe1_e150",
        description="RPE1 E150 screen, gem-group halves",
        split_axis="gem-group split",
        modality="CRISPRi",
        source="scPerturb Zenodo 7041849",
    ),
    "tahoe100m": HalfDataset(
        name="tahoe100m",
        description="Tahoe-100M drug x cell-line context, plate3 "
                    "sample-level halves",
        split_axis="sample split",
        modality="chemical perturbation",
        source="Tahoe-100M",
    ),
}


def load_halves(
    path: str | Path,
    min_cells: int = 25,
    drop_non_targeting: bool = True,
    normalize_by_n: bool | None = None,
):
    """Load a frozen halves npz.

    Returns (half_a, half_b, perts, genes) with shapes
    (n_pert, n_genes), (n_pert, n_genes), (n_pert,), (n_genes,).

    Real frozen archives store summed expression in ``S`` and cell counts in
    ``N``; when ``normalize_by_n`` is True (or auto-detected from un-divided
    sums in ``S`` when ``N`` is present), perturbations are filtered to
    ``N.min(axis=0) >= min_cells`` (and optionally non-targeting excluded)
    and ``S`` is divided by ``N``.
    """
    z = np.load(path, allow_pickle=True)
    s = np.asarray(z["S"], dtype=np.float64)
    if s.shape[0] != 2:
        raise ValueError(f"expected S with leading axis 2, got {s.shape}")
    perts = np.asarray(z["perts"])
    genes = np.asarray(z["genes"])
    if "N" in z.files:
        n = np.asarray(z["N"], dtype=np.float64)
        should_norm = (
            normalize_by_n
            if normalize_by_n is not None
            else (np.abs(s).max() > 50.0 and n.max() > 0)
        )
        if should_norm:
            keep = n.min(axis=0) >= min_cells
            if drop_non_targeting:
                keep = keep & (perts.astype(str) != "non-targeting")
            if not np.all(keep):
                s = s[:, keep]
                n = n[:, keep]
                perts = perts[keep]
            s /= n[:, :, None]
    return s[0], s[1], perts, genes


def list_datasets() -> list[HalfDataset]:
    return list(DATASETS.values())
