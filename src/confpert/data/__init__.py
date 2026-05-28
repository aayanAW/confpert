"""Data loaders for Perturb-seq datasets used in HetPert.

Norman 2019 (K562, single+double perturbations) is the headline diagnostic substrate;
Replogle 2022 K562 essential is the second; Adamson 2016 third; Tahoe-100M optional.

Each loader returns a `PerturbDataset` with:
  X_ctrl: [n_ctrl_cells, d_genes]      -- control population
  X_pert: dict[pert_name -> [n_p_cells, d_genes]]  -- per-perturbation populations
  gene_names: list[str]                 -- column names
  metadata: dict                        -- source dataset name, n_top_genes, etc.
"""
from .adamson import load_adamson
from .norman import NormanDataset, PerturbDataset, load_norman
from .replogle import load_replogle
from .replogle_rpe1 import load_replogle_rpe1
from .tahoe import load_tahoe
from ._phase2_stubs import (
    load_frangieh, load_schmidt, load_datlinger,
    load_mcfaline_figueroa, load_walker, load_lara_astiaso,
)

__all__ = ["NormanDataset", "PerturbDataset", "load_norman", "load_replogle",
           "load_adamson", "load_replogle_rpe1", "load_tahoe",
           "load_frangieh", "load_schmidt", "load_datlinger",
           "load_mcfaline_figueroa", "load_walker", "load_lara_astiaso"]
