"""Phase 2 predictor wrapper stubs.

Each new Phase 2 predictor (scGPT, scFoundation, Geneformer, CellFM, GET)
declared here with the universal `predict_samples` interface. The wrapping
logic itself runs Modal-resident in `scripts/modal_launch.py` (added in
Phase 2B alongside actual training entrypoints).

Each stub declares:
  - id (matches preregistration_v2.yaml `k1_v2.predictors[].id`)
  - family (matches K2 v2 H3 family assignment)
  - param_count (for K2 v2 H1 / H2 capacity covariate)
  - corpus_diversity_index (for K2 v2 H2b data-corpus covariate)
  - license (for datasheet)
  - upstream_repo / model_hub_path
  - heavyweight flag (for K2 v2 §1.6 dataset restriction)
  - status: "stub" (Phase 2A) | "wrapped" (Phase 2B) | "results-merged" (Phase 2C+)

This module is import-safe: no heavy deps (transformers, scGPT, etc.) are
imported here. Heavy imports happen inside the Modal function only.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass(frozen=True)
class Phase2PredictorSpec:
    """Static metadata for a Phase 2 predictor.

    `predict_samples` interface (implemented Modal-side in Phase 2B):
        def predict_samples(
            X_ctrl_test: np.ndarray,  # (n_cells, n_genes)
            perturbation: str,
            n_cells: int,
        ) -> np.ndarray:               # (n_cells, n_genes)
    """

    id: str
    family: str
    param_count: int
    corpus_diversity_index: int
    license: str
    upstream_repo: str
    model_hub_path: str
    heavyweight: bool
    status: str = "stub"
    notes: str = ""


SCGPT = Phase2PredictorSpec(
    id="scgpt",
    family="F3_transformer",
    param_count=int(1.0e8),
    corpus_diversity_index=33,
    license="MIT",
    upstream_repo="https://github.com/bowang-lab/scGPT",
    model_hub_path="bowang-lab/scgpt-pretrained",
    heavyweight=True,
    notes=(
        "Bin-tokenised transformer. Phase 2B: implement Modal entrypoint following "
        "the GEARS pattern in modal_launch.py. Requires HF auth for some weights. "
        "Tokenizer + binning step is the trickiest integration point. "
        "Per modal/scgpt-discussion #142 (community), pin transformers<4.46."
    ),
)


SCFOUNDATION = Phase2PredictorSpec(
    id="scfoundation",
    family="F3_transformer",
    param_count=int(1.0e8),
    corpus_diversity_index=50,
    license="Apache-2.0",
    upstream_repo="https://github.com/biomap-research/scFoundation",
    model_hub_path="biomap-research/scfoundation",
    heavyweight=True,
    notes=(
        "Asymmetric transformer-like architecture. Phase 2B: mirror scGPT integration. "
        "xTrimo serving format used internally; HF mirror checkpoints recommended."
    ),
)


GENEFORMER = Phase2PredictorSpec(
    id="geneformer",
    family="F3_transformer",
    param_count=int(1.0e7),
    corpus_diversity_index=30,
    license="Apache-2.0",
    upstream_repo="https://huggingface.co/ctheodoris/Geneformer",
    model_hub_path="ctheodoris/Geneformer",
    heavyweight=True,
    notes=(
        "Rank-value encoding (distinct from scGPT/scFoundation binning). Tests H3 family "
        "ordering: if H3 lands, this distinguishes 'rank-encoding' from 'binning' within "
        "the F3_transformer family. Phase 2B."
    ),
)


CELLFM = Phase2PredictorSpec(
    id="cellfm",
    family="F3_transformer",
    param_count=int(8.0e8),
    corpus_diversity_index=100,
    license="Apache-2.0",
    upstream_repo="https://github.com/biomed-AI/CellFM",
    model_hub_path="biomed-AI/cellfm-800M",
    heavyweight=True,
    notes=(
        "RetNet architecture on MindSpore. 800M params — largest open scFM. Stress-tests "
        "K2 v2 H1 capacity at 4x STATE. Pinning MindSpore + Linux + CUDA combination is "
        "the integration risk. Phase 2B may need to skip if MindSpore install fails on "
        "Modal A100 image. Backup: HF-mirrored PyTorch port if available by then."
    ),
)


GET = Phase2PredictorSpec(
    id="get",
    family="EXCLUDED_H3",
    param_count=int(5.0e7),
    corpus_diversity_index=213,
    license="MIT",
    upstream_repo="https://github.com/GET-Foundation/get_model",
    model_hub_path="get-foundation/get-pretrained",
    heavyweight=True,
    status="deferred_optional",
    notes=(
        "Chromatin-conditioned. Distinct input modality (ATAC + sequence, not RNA only). "
        "Demoted to optional extension per PHASE_2_PLAN.md §14.5 (Modal budget cap). "
        "Run only if budget allows after the must-have 13 predictors complete."
    ),
)


PHASE_2_PREDICTORS: dict[str, Phase2PredictorSpec] = {
    p.id: p for p in [SCGPT, SCFOUNDATION, GENEFORMER, CELLFM, GET]
}


def list_phase2_predictors(include_deferred: bool = True) -> list[Phase2PredictorSpec]:
    """List Phase 2 predictor specs in priority order."""
    out = [SCGPT, SCFOUNDATION, GENEFORMER, CELLFM]
    if include_deferred:
        out.append(GET)
    return out


def heavyweight_dataset_allowlist() -> list[str]:
    """Datasets that heavyweight Phase 2 predictors are allowed to run on,
    per preregistration_v2.yaml §1.6 (k1_v2.heavyweight_restriction).
    """
    return ["norman", "replogle_k562", "replogle_rpe1", "tahoe", "frangieh", "schmidt"]


def is_heavyweight_predictor(predictor_id: str) -> bool:
    """Return True if predictor is in the heavyweight class (subject to restriction)."""
    if predictor_id in {"scgpt", "scfoundation", "geneformer", "state", "cellfm"}:
        return True
    return False
