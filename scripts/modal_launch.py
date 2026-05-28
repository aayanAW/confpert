"""Modal --detach launchers for the heavy ConfPert predictors.

Per PHASE_1_PLAN.md and compute_estimate.md:
  - GEARS uncertainty mode: run cell-gears v0.1.2 with predict() exposing sigma_hat^2.
  - CPA: cpa-tools, retrain per dataset.
  - biolord: latent-optimization predictor, retrain per dataset.
  - sVAE+ / SAMS-VAE: insitro/sams-vae, retrain per dataset.
  - STATE: ArcInstitute/state with SE-600M public checkpoint, inference + calibration only.

All entrypoints use the shared `causeflow-artifacts` Modal volume so we can pull
locally afterward. Detached mode is the default per the user's overnight directive.

Launch:
  modal run --detach scripts/modal_launch.py::gears_calibrate_norman
  modal run --detach scripts/modal_launch.py::cpa_calibrate_norman
"""
from __future__ import annotations

import os
from pathlib import Path

import modal

APP_NAME = "confpert"

# Shared image with torch + scientific stack. Predictor-specific images extend this.
# Modal requires add_local_* to be LAST so image rebuilds don't trigger on local file edits.
_BASE_DEPS = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.4.1",
        "numpy>=1.26",
        "scipy>=1.11",
        "scikit-learn>=1.4",
        "pandas>=2.1",
        "anndata>=0.10",
        "scanpy>=1.10",
        "tqdm>=4.66",
    )
)


def _add_local(img: modal.Image) -> modal.Image:
    return (
        img.add_local_dir(str(Path(__file__).resolve().parent.parent / "src"),
                          remote_path="/root/src")
           .add_local_dir(str(Path(__file__).resolve().parent.parent / "scripts"),
                          remote_path="/root/scripts")
    )


BASE_IMAGE = _add_local(_BASE_DEPS)

GEARS_IMAGE = _add_local(
    _BASE_DEPS.pip_install(
        "torch_geometric==2.7.0",
        "cell-gears==0.1.2",
    )
)

# scvi-tools 1.4 changed callback contracts in ways that make biolord 0.0.3
# return loss=None (TypeError in on_train_batch_end isnan check) and removes
# the scvi._compat module that scgen 2.1.0 imports. Pin scvi-tools<1.4 across
# all three predictors. CPA also needs libtk8.6 for matplotlib's tkinter import
# at cpa-tools import time; MPLBACKEND=Agg avoids the GUI backend negotiation.
#
# CPA: cpa-tools 0.8.8 (latest) requires Python >=3.9,<3.11. Building on the
# default 3.11 base causes pip to resolve to the empty placeholder 0.0.0,
# which is missing rdkit. Use a separate 3.10 base for CPA only.
_CPA_BASE = (
    modal.Image.debian_slim(python_version="3.10")
    .apt_install("libtk8.6")
    .env({"MPLBACKEND": "Agg"})
    .pip_install(
        "torch==2.4.1",
        "numpy>=1.26",
        "scipy>=1.11",
        "scikit-learn>=1.4",
        "pandas>=2.1",
        "anndata>=0.10",
        "scanpy>=1.10",
        "tqdm>=4.66",
        "rdkit",
    )
)
CPA_IMAGE = _add_local(
    _CPA_BASE.pip_install(
        "cpa-tools==0.8.8",
        "scvi-tools>=0.20.3,<1.0",
        # pyarrow >=14 removed PyExtensionType which scvi-tools 0.20 imports.
        "pyarrow<14",
    )
)

# scgen 2.1.0 on PyPI imports `scvi._compat` and `LossRecorder`, both removed
# from scvi-tools 1.3+. scgen master (commit d79e1f04, Oct 2023) was patched to
# use `from typing import Literal` and `LossOutput` instead — install from git.
# Modal's debian_slim image lacks git, so apt_install it before pip-from-git.
SCGEN_IMAGE = _add_local(
    _BASE_DEPS.apt_install("git")
              .env({"MPLBACKEND": "Agg"})
              .pip_install(
                  "git+https://github.com/theislab/scgen.git",
                  "scvi-tools>=1.2,<1.4",
              )
)

# biolord 0.0.3 on PyPI returns loss=None on scvi-tools 1.3+. biolord master is
# at the same version 0.0.3 but with later patches; install from git for compat.
BIOLORD_IMAGE = _add_local(
    _BASE_DEPS.apt_install("git")
              .env({"MPLBACKEND": "Agg"})
              .pip_install(
                  "git+https://github.com/nitzanlab/biolord.git",
                  "scvi-tools>=1.2,<1.4",
              )
)

# sVAE+ / SAMS-VAE: insitro/sams-vae has no PyPI and a broken setup.cfg
# (declares `where = src` but package is at top level). pip install from git
# installs zero packages. Workaround: git clone in the image and add the
# top-level repo path to PYTHONPATH at function time.
SVAEPLUS_IMAGE = _add_local(
    _BASE_DEPS.apt_install("git")
              .env({"MPLBACKEND": "Agg", "PYTHONPATH": "/root/src:/opt/sams-vae"})
              .pip_install(
                  "pyro-ppl>=1.8",
                  "pytorch-lightning>=2.0",
                  "tqdm",
              )
              .run_commands(
                  "git clone --depth 1 https://github.com/insitro/sams-vae.git /opt/sams-vae"
              )
)

app = modal.App(APP_NAME, image=BASE_IMAGE)
volume = modal.Volume.from_name("causeflow-artifacts", create_if_missing=True)
VOLUME_MOUNT = "/artifacts"


def _setup_env():
    import sys
    sys.path.insert(0, "/root/src")
    os.environ.setdefault("PYTHONPATH", "/root/src")


# ---------------------------------------------------------------------------
# GEARS in uncertainty mode (the actually-distributional predict() path)
# ---------------------------------------------------------------------------


@app.function(
    gpu="A100",
    timeout=60 * 60 * 2,
    volumes={VOLUME_MOUNT: volume},
    image=GEARS_IMAGE,
)
def _gears_calibrate_fn(dataset: str = "norman", epochs: int = 5,
                        hidden_size: int = 64, n_eval_perts: int = 30):
    print(f"[gears confpert] FUNCTION ENTRY dataset={dataset} epochs={epochs}",
          flush=True)
    """Train GEARS with uncertainty head exposed; calibrate per-discrepancy
    against held-out perturbations on the chosen dataset.
    """
    _setup_env()
    import json
    import time
    import warnings
    warnings.filterwarnings("ignore")

    import numpy as np
    import torch
    import sys
    sys.path.insert(0, "/root/src")

    from gears import PertData, GEARS
    from confpert.metrics import SCORES, evaluate_six
    from confpert.conformal import PerturbationConformal

    print(f"[gears confpert] CUDA={torch.cuda.is_available()} dataset={dataset}",
          flush=True)

    # Pandas Series.nonzero() patch for newer pandas (HetPert known issue)
    import pandas as _pd
    if not hasattr(_pd.Series, "_nonzero_patched"):
        _pd.Series.nonzero = lambda self: self.to_numpy().nonzero()
        _pd.Series._nonzero_patched = True

    gears_data_dir = "/root/gears_data"
    os.makedirs(gears_data_dir, exist_ok=True)
    pd_obj = PertData(gears_data_dir)
    print(f"[gears confpert] loading {dataset} ...", flush=True)
    pd_obj.load(data_name=dataset)
    pd_obj.prepare_split(split="simulation", seed=1)
    pd_obj.get_dataloader(batch_size=64, test_batch_size=64)

    print(f"[gears confpert] training GEARS uncertainty=True for {epochs} ep ...",
          flush=True)
    gears_model = GEARS(pd_obj, device="cuda")
    gears_model.model_initialize(hidden_size=hidden_size, uncertainty=True)
    t0 = time.time()
    gears_model.train(epochs=epochs, lr=1e-3)
    fit_t = time.time() - t0
    print(f"[gears confpert] trained in {fit_t:.0f}s", flush=True)

    # Eval distributional metrics on top-N test perturbations.
    # Prefer subgroup["test_subgroup"] (Norman simulation-split has it); fall
    # back to set2conditions["test"] (Adamson and other datasets where
    # subgroup is None).
    test_pert_names: list[str] = []
    if pd_obj.subgroup and pd_obj.subgroup.get("test_subgroup"):
        for cat in pd_obj.subgroup["test_subgroup"].values():
            test_pert_names.extend(cat)
    elif hasattr(pd_obj, "set2conditions") and "test" in pd_obj.set2conditions:
        test_pert_names = list(pd_obj.set2conditions["test"])
    print(f"[gears confpert] test pert candidates: {len(test_pert_names)}",
          flush=True)

    pert_to_n = {}
    for p_name in test_pert_names:
        # Exact "ctrl" only — Adamson uses "GENE+ctrl" single-pert format
        # which would be wrongly filtered by `"ctrl" in p_name.lower()`.
        if p_name.lower() in {"ctrl", "control"}:
            continue
        mask = pd_obj.adata.obs["condition"] == p_name
        pert_to_n[p_name] = int(mask.sum())
    top_perts = sorted(pert_to_n.items(), key=lambda x: -x[1])[:n_eval_perts]
    print(f"[gears confpert] eval top-{n_eval_perts} non-ctrl perts", flush=True)

    pred_pops = {}
    obs_pops = {}
    pert_t0 = time.time()
    for pi, (p_name, _) in enumerate(top_perts):
        obs_mask = (pd_obj.adata.obs["condition"] == p_name).values
        X_obs = pd_obj.adata.X[obs_mask]
        if hasattr(X_obs, "toarray"):
            X_obs = X_obs.toarray()
        X_obs = np.asarray(X_obs).astype(np.float32)
        if X_obs.shape[0] < 10:
            continue
        try:
            # Strip "ctrl" from single-perturbation format "GENE+ctrl" (Adamson
            # convention); cell-gears' predict() expects [["GENE"]] for singles
            # and rejects "ctrl" with "ctrl is not in the perturbation graph".
            pert_list = [g for g in p_name.split("+") if g.lower() != "ctrl"]
            if not pert_list:
                continue
            t_pred = time.time()
            print(f"[gears confpert] {pi+1}/{len(top_perts)} predict {p_name} "
                  f"(pert_list={pert_list}) ...", flush=True)
            pred = gears_model.predict([pert_list])
            print(f"[gears confpert] {pi+1}/{len(top_perts)} predict {p_name} "
                  f"done in {time.time()-t_pred:.1f}s", flush=True)
            # Uncertainty-mode predict returns dict with both 'mean' and 'sigma'
            if isinstance(pred, tuple) and len(pred) == 2:
                pred_mean, pred_sigma = pred
            else:
                pred_mean = pred
                pred_sigma = None
            key = p_name if p_name in pred_mean else (
                "_".join(pert_list) if "_".join(pert_list) in pred_mean
                else list(pred_mean.values())[0]
            )
            X_pred_mean = (np.asarray(pred_mean[key]).flatten()
                           if isinstance(pred_mean, dict) else np.asarray(key).flatten()
                           ).astype(np.float32)
            if pred_sigma is not None and isinstance(pred_sigma, dict):
                sig_key = (key if key in pred_sigma
                           else list(pred_sigma.values())[0])
                X_pred_sigma = np.sqrt(np.maximum(
                    np.asarray(pred_sigma[sig_key]).flatten(), 1e-9
                )).astype(np.float32)
            else:
                X_pred_sigma = np.ones_like(X_pred_mean) * 1e-3
            # Sample marginal Gaussian per gene
            rng = np.random.RandomState(42)
            n_cells = X_obs.shape[0]
            eps = rng.randn(n_cells, X_pred_mean.shape[0]).astype(np.float32)
            X_pred_samples = (X_pred_mean[None, :]
                              + X_pred_sigma[None, :] * eps).astype(np.float32)
            if X_pred_samples.shape[1] != X_obs.shape[1]:
                continue
            pred_pops[p_name] = X_pred_samples
            obs_pops[p_name] = X_obs
        except Exception as e:
            print(f"  predict failed for {p_name}: {type(e).__name__}: {e}", flush=True)

    # Calibrate each of the six discrepancies via PerturbationConformal
    pert_keys = list(pred_pops.keys())
    rng = np.random.RandomState(0)
    rng.shuffle(pert_keys)
    n_split = max(2, len(pert_keys) // 2)
    calib_keys = pert_keys[:n_split]
    test_keys = pert_keys[n_split:]
    pred_calib = [pred_pops[k] for k in calib_keys]
    obs_calib = [obs_pops[k] for k in calib_keys]
    pred_test_list = [pred_pops[k] for k in test_keys]
    obs_test_list = [obs_pops[k] for k in test_keys]

    out = {"config": {"dataset": dataset, "epochs": epochs, "hidden_size": hidden_size,
                      "n_eval_perts": n_eval_perts, "predictor": "gears_uncertainty"},
           "fit_sec": fit_t,
           "n_pred_pops": len(pred_pops),
           "calibration_results": {}}
    for alpha in [0.05, 0.10, 0.20]:
        out["calibration_results"][f"alpha_{alpha}"] = {}
        for score_name, score_fn in SCORES.items():
            try:
                pc = PerturbationConformal(score_fn=score_fn, alpha=alpha)
                pc.calibrate(pred_calib, obs_calib)
                cov = pc.coverage(pred_test_list, obs_test_list)
                out["calibration_results"][f"alpha_{alpha}"][score_name] = cov
            except Exception as e:
                out["calibration_results"][f"alpha_{alpha}"][score_name] = {
                    "error": f"{type(e).__name__}: {e}"
                }

    out_path = f"{VOLUME_MOUNT}/confpert_gears_{dataset}.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=float)
    volume.commit()
    print(f"[gears confpert] -> {out_path}", flush=True)
    return out


@app.local_entrypoint()
def gears_calibrate_norman():
    """GEARS uncertainty + ConfPert calibration on Norman. ~$5, ~30 min."""
    _gears_calibrate_fn.remote(dataset="norman", epochs=5)


@app.local_entrypoint()
def gears_calibrate(dataset: str = "norman", epochs: int = 5):
    """Parameterized GEARS entrypoint. dataset must match cell-gears built-in
    (norman, adamson, dixit). Replogle K562/RPE1 are NOT in cell-gears' built-in
    download list and would require a custom PertData prep.
    """
    _gears_calibrate_fn.remote(dataset=dataset, epochs=epochs)


# ---------------------------------------------------------------------------
# GEARS-uncertainty on Replogle K562/RPE1 via custom PertData.new_data_process
# ---------------------------------------------------------------------------


@app.function(
    gpu="A100",
    timeout=60 * 60 * 3,
    volumes={VOLUME_MOUNT: volume},
    image=GEARS_IMAGE,
)
def _gears_calibrate_custom_fn(dataset: str = "replogle_k562", epochs: int = 5,
                                n_eval_perts: int = 30,
                                n_top_perturbations: int = 50):
    """Train GEARS on a non-built-in dataset (Replogle K562/RPE1) by building
    a custom PertData via new_data_process. Maps our on-volume Replogle h5ad
    obs schema to cell-gears' expected schema (`condition` with 'ctrl' / 'GENE'
    / 'GENE1+GENE2', `cell_type`, var.gene_name).
    """
    print(f"[gears-custom confpert] FUNCTION ENTRY dataset={dataset} epochs={epochs}",
          flush=True)
    _setup_env()
    import os
    import time
    import warnings
    warnings.filterwarnings("ignore")

    import numpy as np
    import pandas as pd
    import scanpy as sc
    import torch
    import sys
    sys.path.insert(0, "/root/src")

    from gears import PertData, GEARS
    from confpert.metrics import SCORES, evaluate_six
    from confpert.conformal import PerturbationConformal

    print(f"[gears-custom confpert] CUDA={torch.cuda.is_available()} dataset={dataset}",
          flush=True)

    # Pandas Series.nonzero() patch for newer pandas
    import pandas as _pd
    if not hasattr(_pd.Series, "_nonzero_patched"):
        _pd.Series.nonzero = lambda self: self.to_numpy().nonzero()
        _pd.Series._nonzero_patched = True

    # Load on-volume h5ad and adapt schema to cell-gears expectations
    h5ad_path = _h5ad_path_on_volume(dataset)
    adata = sc.read_h5ad(h5ad_path)
    print(f"[gears-custom] raw shape: {adata.shape}; obs cols: "
          f"{list(adata.obs.columns)[:8]}...", flush=True)

    # Map perturbation column → 'condition' (cell-gears expects 'GENE' for
    # singles or 'ctrl' for control). Replogle uses 'control' / 'GENE'.
    pert_col = None
    for c in ("perturbation", "gene", "gene_id"):
        if c in adata.obs.columns:
            pert_col = c
            break
    if pert_col is None:
        raise ValueError(f"No perturbation column in {h5ad_path}")
    labels = adata.obs[pert_col].astype(str)
    # Normalize: control → 'ctrl', single perts → 'GENE+ctrl' (cell-gears
    # convention for singles in datasets that include doubles).
    is_ctrl = labels.str.lower().isin({"control", "ctrl", "non-targeting"})
    new_cond = labels.where(~is_ctrl, "ctrl")
    adata.obs["condition"] = new_cond.where(is_ctrl, new_cond + "+ctrl")
    adata.obs["cell_type"] = "K562" if "k562" in dataset.lower() else "RPE1"

    # Subsample to top n_top_perturbations + control to keep DE + PyG build tractable
    counts = adata.obs[pert_col][~is_ctrl].value_counts()
    top_perts = counts.head(n_top_perturbations).index.tolist()
    keep_mask = is_ctrl | adata.obs[pert_col].isin(top_perts)
    adata = adata[keep_mask.values].copy()
    print(f"[gears-custom] after subset: shape={adata.shape}, n_perts="
          f"{adata.obs['condition'].nunique()}", flush=True)

    # gene_name in var for cell-gears
    if "gene_name" not in adata.var.columns:
        adata.var["gene_name"] = adata.var.index.astype(str)

    # cell-gears expects HVG selection + log1p before training is helpful
    if "log1p" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    if adata.shape[1] > 5000:
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
    # cell-gears' get_dropout_non_zero_genes calls X.toarray() so adata.X must
    # be sparse going into new_data_process. After highly_variable_genes(subset=True)
    # the X is densified; convert back.
    from scipy.sparse import issparse, csr_matrix
    if not issparse(adata.X):
        adata.X = csr_matrix(adata.X)
    print(f"[gears-custom] post-preproc: {adata.shape}, X sparse={issparse(adata.X)}", flush=True)

    # Build cell-gears PertData via new_data_process
    gears_data_dir = "/root/gears_data_custom"
    os.makedirs(gears_data_dir, exist_ok=True)
    pd_obj = PertData(gears_data_dir)
    print(f"[gears-custom] new_data_process for {dataset} ...", flush=True)
    t_prep = time.time()
    pd_obj.new_data_process(dataset_name=dataset, adata=adata,
                             skip_calc_de=False)
    print(f"[gears-custom] new_data_process done in {time.time()-t_prep:.0f}s",
          flush=True)

    pd_obj.prepare_split(split="simulation", seed=1)
    pd_obj.get_dataloader(batch_size=64, test_batch_size=64)

    print(f"[gears-custom] training GEARS uncertainty=True for {epochs} ep ...",
          flush=True)
    gears_model = GEARS(pd_obj, device="cuda")
    gears_model.model_initialize(hidden_size=64, uncertainty=True)
    t0 = time.time()
    gears_model.train(epochs=epochs, lr=1e-3)
    fit_t = time.time() - t0
    print(f"[gears-custom] trained in {fit_t:.0f}s", flush=True)

    # Eval (same code path as built-in version)
    test_pert_names: list[str] = []
    if pd_obj.subgroup and pd_obj.subgroup.get("test_subgroup"):
        for cat in pd_obj.subgroup["test_subgroup"].values():
            test_pert_names.extend(cat)
    elif hasattr(pd_obj, "set2conditions") and "test" in pd_obj.set2conditions:
        test_pert_names = list(pd_obj.set2conditions["test"])
    print(f"[gears-custom] test pert candidates: {len(test_pert_names)}",
          flush=True)

    pert_to_n = {}
    for pname in test_pert_names:
        if pname.lower() in {"ctrl", "control"}:
            continue
        mask = pd_obj.adata.obs["condition"] == pname
        pert_to_n[pname] = int(mask.sum())
    top_test = sorted(pert_to_n.items(), key=lambda x: -x[1])[:n_eval_perts]
    print(f"[gears-custom] eval top-{n_eval_perts} non-ctrl perts", flush=True)

    pred_pops, obs_pops = {}, {}
    for pi, (pname, _) in enumerate(top_test):
        obs_mask = (pd_obj.adata.obs["condition"] == pname).values
        X_obs = pd_obj.adata.X[obs_mask]
        if hasattr(X_obs, "toarray"):
            X_obs = X_obs.toarray()
        X_obs = np.asarray(X_obs).astype(np.float32)
        if X_obs.shape[0] < 10:
            continue
        try:
            pert_genes = [g for g in pname.split("+") if g.lower() != "ctrl"]
            if not pert_genes:
                continue
            t_pred = time.time()
            pred = gears_model.predict([pert_genes])
            print(f"[gears-custom] {pi+1}/{len(top_test)} {pname}: predicted "
                  f"in {time.time()-t_pred:.1f}s", flush=True)
            if isinstance(pred, tuple) and len(pred) == 2:
                pred_mean, pred_sigma = pred
            else:
                pred_mean = pred
                pred_sigma = None
            key = pname if isinstance(pred_mean, dict) and pname in pred_mean else (
                "_".join(pert_genes)
                if isinstance(pred_mean, dict) and "_".join(pert_genes) in pred_mean
                else (list(pred_mean.values())[0] if isinstance(pred_mean, dict) else pred_mean)
            )
            X_pred_mean = (np.asarray(pred_mean[key]).flatten()
                           if isinstance(pred_mean, dict) and key in pred_mean
                           else np.asarray(key).flatten()
                           ).astype(np.float32)
            if pred_sigma is not None and isinstance(pred_sigma, dict):
                sig_key = (key if key in pred_sigma
                           else list(pred_sigma.values())[0])
                X_pred_sigma = np.sqrt(np.maximum(
                    np.asarray(pred_sigma[sig_key]).flatten(), 1e-9
                )).astype(np.float32)
            else:
                X_pred_sigma = np.ones_like(X_pred_mean) * 1e-3
            rng = np.random.RandomState(42)
            n_cells = X_obs.shape[0]
            eps = rng.randn(n_cells, X_pred_mean.shape[0]).astype(np.float32)
            X_pred_samples = (X_pred_mean[None, :]
                              + X_pred_sigma[None, :] * eps).astype(np.float32)
            if X_pred_samples.shape[1] != X_obs.shape[1]:
                continue
            pred_pops[pname] = X_pred_samples
            obs_pops[pname] = X_obs
        except Exception as e:
            print(f"  [gears-custom] predict failed for {pname}: "
                  f"{type(e).__name__}: {e}", flush=True)

    print(f"[gears-custom] paired {len(pred_pops)} perts; running conformal ...",
          flush=True)
    out = _run_conformal_eval(pred_pops, obs_pops, dataset, "gears_uncertainty",
                               fit_t, extra_config={
                                   "epochs": epochs,
                                   "n_eval_perts": n_eval_perts,
                                   "custom_pert_data": True})
    volume.commit()
    print(f"[gears-custom] -> {VOLUME_MOUNT}/confpert_gears_uncertainty_{dataset}.json",
          flush=True)
    return out


@app.local_entrypoint()
def gears_calibrate_custom(dataset: str = "replogle_k562", epochs: int = 5):
    """GEARS on Replogle K562/RPE1 via custom PertData prep. Slower than
    built-in datasets due to new_data_process (DE + PyG graph build)."""
    _gears_calibrate_custom_fn.remote(dataset=dataset, epochs=epochs)


# ---------------------------------------------------------------------------
# Shared helpers for sample-producing predictors
# ---------------------------------------------------------------------------


def _h5ad_path_on_volume(dataset: str) -> str:
    """Return the absolute path of the dataset's h5ad on the mounted volume.

    Note: the on-volume Norman / Replogle K562 h5ads have hetpert_-free names
    in /artifacts/data/ (verified via `modal volume ls`).
    """
    name = {
        "norman": "norman_2019.h5ad",
        "replogle_k562": "replogle_k562_essential.h5ad",
        "replogle_rpe1": "replogle_rpe1.h5ad",
        "adamson": "adamson_2016.h5ad",
        "tahoe": "tahoe_subset.h5ad",
        # Phase 2C datasets (session 6 prestaged via modal_launch phase2_dataset_prestage)
        "frangieh": "frangieh.h5ad",
        "datlinger": "datlinger.h5ad",
        "schmidt": "schmidt.h5ad",
        "mcfaline_figueroa": "mcfaline_figueroa.h5ad",
    }.get(dataset)
    if name is None:
        raise KeyError(
            f"No on-volume h5ad for {dataset}; expected one of "
            f"norman / replogle_k562 / replogle_rpe1 / adamson / tahoe / "
            f"frangieh / datlinger / schmidt / mcfaline_figueroa."
        )
    return f"{VOLUME_MOUNT}/data/{name}"


# Tahoe uses `drug` as the perturbation column and `DMSO` as the control marker
# (vs Norman/Replogle/Adamson which use a gene-name column with regexable
# control patterns). Branch on dataset name.
TAHOE_CTRL_REGEX = "dmso|vehicle|^ctrl$"


def _find_ctrl_value(labels) -> str | None:
    """Return the in-data string used for control / vehicle perturbation, or None.

    Recognises:
      - 'ctrl', 'control'                (Norman, Replogle, Adamson)
      - 'non-targeting', 'NTC',
        'non_target', 'nontarget',
        'scramble', 'safe'               (Replogle and Replogle-pattern guides)
      - 'DMSO', 'DMSO_TF'                (Tahoe-100M)
      - 'vehicle'                        (general drug-screen convention)
    """
    for v in labels.unique():
        vl = str(v).lower()
        if "ctrl" in vl or vl == "control":
            return v
        if "non-target" in vl or "non_target" in vl or "nontarget" in vl:
            return v
        if vl == "ntc" or vl == "scramble" or vl == "safe":
            return v
        if vl.startswith("dmso") or vl == "vehicle":
            return v
    return None


def _load_perturb_ds(dataset: str, n_top_perturbations: int = 50, n_hvg: int = 512,
                     min_cells_per_pert: int = 30):
    """Load a Norman-style PerturbDataset on Modal using the on-volume h5ad.

    Mirrors confpert.data.load_norman preprocessing but reads from /artifacts/data/.
    Returns a `PerturbDataset` plus the underlying scanpy AnnData (post-HVG-subset)
    for predictors (CPA, scGen, biolord) that want to reuse the AnnData object.
    """
    import scanpy as sc
    from confpert.data.norman import CTRL_REGEX, PerturbDataset

    h5ad = _h5ad_path_on_volume(dataset)
    adata = sc.read_h5ad(h5ad)

    if dataset in {"tahoe", "mcfaline_figueroa"}:
        # Chemical-perturbation schema: drug column, DMSO/vehicle controls.
        col = None
        for c in ("drug", "compound", "treatment", "feature", "drug_name"):
            if c in adata.obs.columns:
                col = c
                break
        ctrl_regex = TAHOE_CTRL_REGEX
    else:
        col = None
        for c in ("perturbation", "guide_id", "gene_id", "gene", "guide_identity",
                  "perturbation_name", "condition", "gene_name"):
            if c in adata.obs.columns:
                col = c
                break
        ctrl_regex = CTRL_REGEX
    if col is None:
        raise ValueError(f"No perturbation column in {h5ad}: cols={list(adata.obs.columns)}")

    labels = adata.obs[col].astype(str)
    is_ctrl = labels.str.lower().str.contains(ctrl_regex, regex=True)
    pert_labels = labels[~is_ctrl]
    counts = pert_labels.value_counts()
    counts = counts[counts >= min_cells_per_pert]
    top_perts = counts.head(n_top_perturbations).index.tolist()

    keep_mask = is_ctrl | labels.isin(top_perts)
    adata = adata[keep_mask.values].copy()

    sc.pp.filter_cells(adata, min_counts=200)
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    adata = adata[:, adata.var["highly_variable"]].copy()

    import numpy as np
    X_full = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    X_full = X_full.astype(np.float32)

    labels_post = adata.obs[col].astype(str)
    is_ctrl_post = labels_post.str.lower().str.contains(ctrl_regex, regex=True)
    X_ctrl = X_full[is_ctrl_post.values]

    X_pert = {}
    for p in top_perts:
        mask = (labels_post == p).values
        if mask.sum() >= min_cells_per_pert:
            X_pert[p] = X_full[mask]

    ds = PerturbDataset(
        X_ctrl=X_ctrl, X_pert=X_pert,
        gene_names=adata.var_names.tolist(),
        metadata={
            "source": f"on-volume {dataset}",
            "perturbation_column": col,
            "n_top_perturbations_returned": len(X_pert),
            "n_ctrl": int(X_ctrl.shape[0]),
            "n_hvg": X_full.shape[1],
        })
    return ds, adata, col, is_ctrl_post.values


def _run_conformal_eval(pred_pops: dict, obs_pops: dict, dataset: str,
                        predictor_name: str, fit_sec: float,
                        extra_config: dict | None = None) -> dict:
    """Run PerturbationConformal calibration on six discrepancies x three alphas.

    Splits the perturbations 50/50 into calib/test. Returns the JSON dict in the
    same shape as `_gears_calibrate_fn`'s output, ready for merge_modal_results.py.
    """
    import json
    import numpy as np
    from confpert.metrics import SCORES
    from confpert.conformal import PerturbationConformal

    pert_keys = list(pred_pops.keys())
    rng = np.random.RandomState(0)
    rng.shuffle(pert_keys)
    n_split = max(2, len(pert_keys) // 2)
    calib_keys = pert_keys[:n_split]
    test_keys = pert_keys[n_split:]
    pred_calib = [pred_pops[k] for k in calib_keys]
    obs_calib = [obs_pops[k] for k in calib_keys]
    pred_test_list = [pred_pops[k] for k in test_keys]
    obs_test_list = [obs_pops[k] for k in test_keys]

    config = {"dataset": dataset, "predictor": predictor_name}
    if extra_config:
        config.update(extra_config)

    out = {"config": config, "fit_sec": float(fit_sec),
           "n_pred_pops": int(len(pred_pops)),
           "calibration_results": {}}
    for alpha in [0.05, 0.10, 0.20]:
        out["calibration_results"][f"alpha_{alpha}"] = {}
        for score_name, score_fn in SCORES.items():
            try:
                pc = PerturbationConformal(score_fn=score_fn, alpha=alpha)
                pc.calibrate(pred_calib, obs_calib)
                cov = pc.coverage(pred_test_list, obs_test_list)
                out["calibration_results"][f"alpha_{alpha}"][score_name] = cov
            except Exception as e:
                out["calibration_results"][f"alpha_{alpha}"][score_name] = {
                    "error": f"{type(e).__name__}: {e}"
                }

    out_path = f"{VOLUME_MOUNT}/confpert_{predictor_name}_{dataset}.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=float)
    return out


def _load_cross_cell_line_data(n_top_per_dataset: int = 200,
                                 train_frac: float = 0.5,
                                 calib_frac: float = 0.25,
                                 seed: int = 42):
    """K562 -> RPE1 cross-cell-line split for heavyweight Modal predictors.

    Mirrors scripts/run_k1_cross_cell_line.py's split logic but operates on
    the on-volume Replogle h5ads. Returns a K562-only AnnData ready for
    setup_anndata / .train(adata) plus held-out arrays for both the calib
    (K562) and test (RPE1) arms in the gene-intersected space.
    """
    import scanpy as sc
    import numpy as np
    import sys
    sys.path.insert(0, "/root/src")
    from confpert.data import load_replogle, load_replogle_rpe1

    print(f"[xcl-load] loading K562 (n_top={n_top_per_dataset}) ...", flush=True)
    ds_k562 = load_replogle(
        h5ad_path="/artifacts/data/replogle_k562_essential.h5ad",
        n_top_perturbations=n_top_per_dataset,
        chunked=True,
    )
    print(f"[xcl-load]   K562: n_ctrl={ds_k562.X_ctrl.shape[0]} "
          f"n_perts={len(ds_k562.X_pert)} d={ds_k562.X_ctrl.shape[1]}", flush=True)

    print(f"[xcl-load] loading RPE1 (n_top={n_top_per_dataset}) ...", flush=True)
    ds_rpe1 = load_replogle_rpe1(
        h5ad_path="/artifacts/data/replogle_rpe1.h5ad",
        n_top_perturbations=n_top_per_dataset,
        chunked=True,
    )
    print(f"[xcl-load]   RPE1: n_ctrl={ds_rpe1.X_ctrl.shape[0]} "
          f"n_perts={len(ds_rpe1.X_pert)} d={ds_rpe1.X_ctrl.shape[1]}", flush=True)

    # Gene intersection
    g_k562 = list(ds_k562.gene_names)
    g_rpe1_set = set(ds_rpe1.gene_names)
    keep_k562 = [i for i, g in enumerate(g_k562) if g in g_rpe1_set]
    if len(keep_k562) < 50:
        raise RuntimeError(f"Only {len(keep_k562)} genes overlap; cross-cell-line needs >=50")
    train_idx = np.asarray(keep_k562, dtype=np.int64)
    train_genes = [g_k562[i] for i in keep_k562]
    test_name_to_idx = {g: i for i, g in enumerate(ds_rpe1.gene_names)}
    test_idx = np.asarray([test_name_to_idx[g] for g in train_genes], dtype=np.int64)
    print(f"[xcl-load]   gene intersection: {len(train_genes)}", flush=True)

    X_ctrl_k562 = ds_k562.X_ctrl[:, train_idx].astype(np.float32)
    X_ctrl_rpe1 = ds_rpe1.X_ctrl[:, test_idx].astype(np.float32)

    # Cell-level split on K562 ctrl cells (50% train, rest calib)
    rng = np.random.RandomState(seed)
    n_ctrl = X_ctrl_k562.shape[0]
    perm = rng.permutation(n_ctrl)
    n_train = int(n_ctrl * train_frac)
    X_ctrl_train = X_ctrl_k562[perm[:n_train]]
    X_ctrl_calib = X_ctrl_k562[perm[n_train:]]

    common_perts = sorted(set(ds_k562.X_pert.keys()) & set(ds_rpe1.X_pert.keys()))
    print(f"[xcl-load]   K562 ∩ RPE1 perts (n_top={n_top_per_dataset}): "
          f"{len(common_perts)}", flush=True)
    if len(common_perts) < 5:
        raise RuntimeError(f"Only {len(common_perts)} pert overlap; need >=5")

    X_pert_train_dict, X_pert_calib_dict, X_pert_test_dict = {}, {}, {}
    for p, X_p in ds_k562.X_pert.items():
        X_p = X_p[:, train_idx].astype(np.float32)
        n_p = X_p.shape[0]
        perm_p = rng.permutation(n_p)
        n_train_p = int(n_p * train_frac)
        n_calib_p = int(n_p * calib_frac)
        X_pert_train_dict[p] = X_p[perm_p[:n_train_p]]
        if p in common_perts:
            X_pert_calib_dict[p] = X_p[perm_p[n_train_p:n_train_p + n_calib_p]]

    for p in common_perts:
        X_p = ds_rpe1.X_pert[p][:, test_idx].astype(np.float32)
        if X_p.shape[0] >= 5:
            X_pert_test_dict[p] = X_p

    common_perts = [p for p in common_perts
                    if p in X_pert_calib_dict and p in X_pert_test_dict
                    and X_pert_calib_dict[p].shape[0] >= 5
                    and X_pert_test_dict[p].shape[0] >= 5]
    print(f"[xcl-load]   common perts with >=5 cells per arm: {len(common_perts)}",
          flush=True)
    if len(common_perts) < 4:
        raise RuntimeError(f"Only {len(common_perts)} usable perts; need >=4")

    # Build K562 train AnnData (ctrl + train perts) for predictor fit
    import anndata as ad
    import pandas as pd
    X_train_list = [X_ctrl_train]
    pert_labels_train = ["non-targeting"] * X_ctrl_train.shape[0]
    for p, X_p in X_pert_train_dict.items():
        X_train_list.append(X_p)
        pert_labels_train.extend([p] * X_p.shape[0])
    X_train_full = np.concatenate(X_train_list, axis=0).astype(np.float32)
    train_adata = ad.AnnData(
        X=X_train_full,
        obs=pd.DataFrame({"perturbation": pert_labels_train,
                          "cell_type": "K562"}),
        var=pd.DataFrame(index=train_genes),
    )
    train_adata.obs_names_make_unique()
    pert_col = "perturbation"
    ctrl_value = "non-targeting"
    print(f"[xcl-load]   train_adata: n={train_adata.shape[0]} d={train_adata.shape[1]} "
          f"n_perts={len(X_pert_train_dict)}", flush=True)

    return (train_adata, pert_col, ctrl_value,
            X_ctrl_calib, X_ctrl_rpe1,
            X_pert_calib_dict, X_pert_test_dict,
            common_perts, train_genes)


def _run_conformal_eval_xcl(pred_calib_pops: dict, obs_calib_pops: dict,
                              pred_test_pops: dict, obs_test_pops: dict,
                              dataset_pair: str, predictor_name: str,
                              fit_sec: float, extra_config: dict | None = None) -> dict:
    """Cross-cell-line conformal eval: K562 calib arm + RPE1 test arm.

    Unlike `_run_conformal_eval`, splits are externally provided (by
    `_load_cross_cell_line_data`) so there is no internal random shuffle.
    """
    import json
    from confpert.metrics import SCORES
    from confpert.conformal import PerturbationConformal

    common = sorted(set(pred_calib_pops.keys()) & set(pred_test_pops.keys())
                    & set(obs_calib_pops.keys()) & set(obs_test_pops.keys()))
    if len(common) < 4:
        raise RuntimeError(f"Only {len(common)} common pert; need >=4 for conformal")
    pred_calib_list = [pred_calib_pops[k] for k in common]
    obs_calib_list = [obs_calib_pops[k] for k in common]
    pred_test_list = [pred_test_pops[k] for k in common]
    obs_test_list = [obs_test_pops[k] for k in common]

    config = {"dataset": dataset_pair, "predictor": predictor_name,
              "split_type": "cross_cell_line"}
    if extra_config:
        config.update(extra_config)

    out = {"config": config, "fit_sec": float(fit_sec),
           "n_pred_pops": len(common),
           "calibration_results": {}}
    for alpha in [0.05, 0.10, 0.20]:
        out["calibration_results"][f"alpha_{alpha}"] = {}
        for score_name, score_fn in SCORES.items():
            try:
                pc = PerturbationConformal(score_fn=score_fn, alpha=alpha)
                pc.calibrate(pred_calib_list, obs_calib_list)
                cov = pc.coverage(pred_test_list, obs_test_list)
                out["calibration_results"][f"alpha_{alpha}"][score_name] = cov
            except Exception as e:
                out["calibration_results"][f"alpha_{alpha}"][score_name] = {
                    "error": f"{type(e).__name__}: {e}"
                }

    out_path = f"{VOLUME_MOUNT}/confpert_{predictor_name}_xcl_{dataset_pair}.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=float)
    return out


# ---------------------------------------------------------------------------
# CPA (Lotfollahi et al. 2023 Mol Syst Biol; cpa-tools)
# Wrapped per lit_notes/cpa_2023.md: predict() returns means; for samples,
# call module.generative(...) to get a Distribution px, then px.sample((S,)).
# ---------------------------------------------------------------------------


@app.function(
    gpu="A100",
    timeout=60 * 60 * 3,
    volumes={VOLUME_MOUNT: volume},
    image=CPA_IMAGE,
)
def _cpa_calibrate_fn(dataset: str = "norman", epochs: int = 50,
                       n_eval_perts: int = 30, n_top_perturbations: int = 50,
                       split_type: str = "within_perturbation"):
    """Train CPA on a Norman-style Perturb-seq dataset, then for each held-out
    perturbation draw S per-cell samples via the generative distribution and
    feed (pred_pop, obs_pop) pairs into ConfPert's PerturbationConformal head.

    Per lit_notes/cpa_2023.md:
      - cpa.CPA.setup_anndata expects perturbation_key, dosage_key, control_group
      - Default predict() returns means via px.loc / px.mu
      - Sample path: model.module.generative(...)["px"].sample((S,))

    split_type='cross_cell_line' trains on Replogle K562 (gene-intersected
    with RPE1) and tests on RPE1 ctrl input, calibrates on K562 held-out.
    """
    print(f"[cpa confpert] FUNCTION ENTRY dataset={dataset} epochs={epochs} "
          f"split={split_type}", flush=True)
    _setup_env()
    import time
    import warnings
    warnings.filterwarnings("ignore")

    import numpy as np
    import torch
    import sys
    sys.path.insert(0, "/root/src")

    import cpa
    from cpa import CPA

    print(f"[cpa confpert] CUDA={torch.cuda.is_available()}", flush=True)

    if split_type == "cross_cell_line":
        (adata, pert_col, ctrl_value, X_ctrl_calib, X_ctrl_rpe1,
         X_pert_calib_dict, X_pert_test_dict, common_perts,
         _gene_names) = _load_cross_cell_line_data()
        print(f"[cpa confpert xcl] train_adata: n={adata.shape[0]}, "
              f"d={adata.shape[1]}, common_perts={len(common_perts)}", flush=True)
        labels = adata.obs[pert_col].astype(str)
        adata.obs["dose_value"] = (labels != ctrl_value).astype(float)
        adata.obs[pert_col] = labels.values
        is_ctrl_mask = (labels == ctrl_value).values
    else:
        ds, adata, pert_col, is_ctrl_mask = _load_perturb_ds(
            dataset, n_top_perturbations=n_top_perturbations
        )
        print(f"[cpa confpert] adata: n={adata.shape[0]}, d={adata.shape[1]}, "
              f"n_perts={len(ds.X_pert)}, n_ctrl={int(is_ctrl_mask.sum())}",
              flush=True)

        # CPA setup: build dose column (1 for any perturbation, 0 for ctrl).
        labels = adata.obs[pert_col].astype(str)
        ctrl_value = _find_ctrl_value(labels)
        if ctrl_value is None:
            raise ValueError(f"No ctrl-like value in {pert_col}: {labels.unique()[:5]}")
        print(f"[cpa confpert] ctrl_value={ctrl_value!r}", flush=True)

        adata.obs["dose_value"] = (labels != ctrl_value).astype(float)
        adata.obs[pert_col] = labels.values

    # cpa-tools ≥0.8 setup_anndata signature
    CPA.setup_anndata(
        adata,
        perturbation_key=pert_col,
        dosage_key="dose_value",
        control_group=ctrl_value,
        batch_key=None,
        is_count_data=False,  # already log1p
    )

    # CPA model with sensible defaults for K562
    model = CPA(
        adata=adata,
        n_latent=64,
        recon_loss="gauss",   # exposes loc + scale at generative time
        doser_type="logsigm",
        n_hidden_encoder=256,
        n_layers_encoder=2,
        n_hidden_decoder=256,
        n_layers_decoder=2,
        use_batch_norm_encoder=True,
        use_layer_norm_encoder=False,
        dropout_rate_encoder=0.1,
        dropout_rate_decoder=0.1,
    )
    print(f"[cpa confpert] training CPA for {epochs} ep ...", flush=True)
    t0 = time.time()
    # cpa-tools 0.8.8's CPATrainingPlan does not accept n_epochs_warmup or
    # n_epochs_pretrain_ae as kwargs (those exist in later cpa versions / scvi
    # 1.x training plans). Use defaults.
    model.train(
        max_epochs=epochs,
        early_stopping_patience=5,
        check_val_every_n_epoch=5,
        batch_size=256,
    )
    fit_t = time.time() - t0
    print(f"[cpa confpert] trained in {fit_t:.0f}s", flush=True)

    # Build per-pert (pred, obs) pairs by drawing samples from the Gaussian head.
    # Skip in cross_cell_line mode (xcl branch below builds its own pred dicts).
    pred_pops, obs_pops = {}, {}
    rng = np.random.RandomState(42)

    if split_type != "cross_cell_line":
        perts_to_eval = sorted(ds.X_pert.keys(),
                               key=lambda p: -ds.X_pert[p].shape[0])[:n_eval_perts]

        # Get a control-cell AnnData slice to feed the encoder (one input cell per
        # output sample; reproduces the per-cell heterogeneity that scGen / CPA
        # rely on for distributional output).
        ctrl_idx = np.where(is_ctrl_mask)[0]

        for p in perts_to_eval:
            n_cells = ds.X_pert[p].shape[0]
            # Sample n_cells control inputs
            if len(ctrl_idx) >= n_cells:
                inp_idx = rng.choice(ctrl_idx, n_cells, replace=False)
            else:
                inp_idx = rng.choice(ctrl_idx, n_cells, replace=True)
            # Build a small AnnData slice with the target perturbation label
            ctrl_slice = adata[inp_idx].copy()
            ctrl_slice.obs[pert_col] = p
            ctrl_slice.obs["dose_value"] = 1.0
            try:
                # Use the high-level model.predict() API. cpa-tools 0.8.8's
                # module.get_expression() (called from predict) does
                # `tensors, _ = self.mixup_data(tensors, alpha=0.0)` which adds the
                # `perts_mixup` field that the lower-level generative() requires.
                # Going through model.predict() handles this; calling generative()
                # directly bypasses mixup_data and crashes with KeyError perts_mixup.
                adata_to_predict = ctrl_slice.copy()
                model.predict(adata=adata_to_predict, n_samples=1, return_mean=False)
                X_pred = adata_to_predict.obsm.get(f"CPA_pred")
                if X_pred is None:
                    # Some cpa-tools versions key by the model class name
                    X_pred = adata_to_predict.obsm.get("CPA_pred")
                if X_pred is None:
                    raise RuntimeError(
                        f"CPA predict left no obsm['CPA_pred']; obsm keys: "
                        f"{list(adata_to_predict.obsm.keys())}"
                    )
                X_pred = np.asarray(X_pred, dtype=np.float32)
                # If 3D (n_samples, n_cells, n_genes) take first sample
                if X_pred.ndim == 3:
                    X_pred = X_pred[0]
                if X_pred.shape[1] != ds.X_pert[p].shape[1]:
                    print(f"  [cpa] dim mismatch for {p}: pred={X_pred.shape}, "
                          f"obs={ds.X_pert[p].shape}", flush=True)
                    continue
                pred_pops[p] = X_pred
                obs_pops[p] = ds.X_pert[p]
            except Exception as e:
                print(f"  [cpa] predict failed for {p}: {type(e).__name__}: {e}",
                      flush=True)

    if split_type == "cross_cell_line":
        # Predict twice per common pert: from K562 ctrl input (calib) +
        # from RPE1 ctrl input (test). Substitute X in a scaffold sliced from
        # train adata so CPA's scvi-tools registry is preserved.
        def _cpa_predict_pop_with_ctrl(X_ctrl_input, p, n_cells, rng_local):
            if X_ctrl_input.shape[0] >= n_cells:
                sel_input = rng_local.choice(X_ctrl_input.shape[0], n_cells,
                                                replace=False)
            else:
                sel_input = rng_local.choice(X_ctrl_input.shape[0], n_cells,
                                                replace=True)
            scaffold_idx = rng_local.choice(adata.shape[0], n_cells,
                                              replace=(adata.shape[0] < n_cells))
            scaffold = adata[scaffold_idx].copy()
            scaffold.X = X_ctrl_input[sel_input].astype(np.float32)
            scaffold.obs[pert_col] = p
            scaffold.obs["dose_value"] = 1.0
            try:
                model.predict(adata=scaffold, n_samples=1, return_mean=False)
            except Exception:
                model.predict(adata=scaffold, n_samples=1)
            X_pred = scaffold.obsm.get("CPA_pred")
            if X_pred is None:
                raise RuntimeError(
                    f"CPA predict left no obsm['CPA_pred']; "
                    f"obsm={list(scaffold.obsm.keys())}"
                )
            X_pred = np.asarray(X_pred, dtype=np.float32)
            if X_pred.ndim == 3:
                X_pred = X_pred[0]
            return X_pred

        pred_calib_pops, obs_calib_pops = {}, {}
        pred_test_pops, obs_test_pops = {}, {}
        for p in common_perts:
            try:
                n_calib = X_pert_calib_dict[p].shape[0]
                n_test = X_pert_test_dict[p].shape[0]
                rng_p = np.random.RandomState(hash(p) % (2**31))
                pred_calib_pops[p] = _cpa_predict_pop_with_ctrl(
                    X_ctrl_calib, p, n_calib, rng_p)
                pred_test_pops[p] = _cpa_predict_pop_with_ctrl(
                    X_ctrl_rpe1, p, n_test, rng_p)
                obs_calib_pops[p] = X_pert_calib_dict[p]
                obs_test_pops[p] = X_pert_test_dict[p]
            except Exception as e:
                print(f"  [cpa xcl] {p}: {type(e).__name__}: {e}", flush=True)
        print(f"[cpa confpert xcl] paired {len(pred_calib_pops)} calib + "
              f"{len(pred_test_pops)} test; running conformal ...", flush=True)
        out = _run_conformal_eval_xcl(
            pred_calib_pops, obs_calib_pops, pred_test_pops, obs_test_pops,
            "replogle_k562_to_rpe1", "cpa", fit_t,
            extra_config={"epochs": epochs}
        )
        volume.commit()
        print(f"[cpa confpert xcl] -> {VOLUME_MOUNT}/confpert_cpa_xcl_"
              f"replogle_k562_to_rpe1.json", flush=True)
        return out

    print(f"[cpa confpert] paired {len(pred_pops)} perts; running conformal ...",
          flush=True)
    out = _run_conformal_eval(pred_pops, obs_pops, dataset, "cpa", fit_t,
                               extra_config={"epochs": epochs,
                                             "n_eval_perts": n_eval_perts})
    volume.commit()
    print(f"[cpa confpert] -> {VOLUME_MOUNT}/confpert_cpa_{dataset}.json",
          flush=True)
    return out


@app.local_entrypoint()
def cpa_calibrate_norman():
    """CPA + ConfPert on Norman. ~$24, ~2-3h on A100."""
    _cpa_calibrate_fn.remote(dataset="norman", epochs=50)


@app.local_entrypoint()
def cpa_calibrate(dataset: str = "norman", epochs: int = 50,
                   split_type: str = "within_perturbation"):
    """Parameterized CPA entrypoint.

    dataset in {norman, replogle_k562, replogle_rpe1, adamson, tahoe} for
    within_perturbation; for split_type=cross_cell_line, dataset is ignored.
    """
    _cpa_calibrate_fn.remote(dataset=dataset, epochs=epochs,
                               split_type=split_type)


# ---------------------------------------------------------------------------
# scGen (Lotfollahi et al. 2019 Nat Methods; theislab/scgen)
# Wrapped per lit_notes/scgen_2019.md: VAE + per-pert latent delta.
# scGen's high-level predict() is celltype-binary; we use the underlying VAE
# directly (encode -> add per-pert delta -> decode) to handle multi-pert Norman.
# ---------------------------------------------------------------------------


@app.function(
    gpu="A100",
    timeout=60 * 60 * 2,
    volumes={VOLUME_MOUNT: volume},
    image=SCGEN_IMAGE,
)
def _scgen_calibrate_fn(dataset: str = "norman", epochs: int = 100,
                         n_eval_perts: int = 30, n_top_perturbations: int = 50,
                         split_type: str = "within_perturbation"):
    """Train scGen on a Norman-style dataset, compute per-pert latent deltas,
    and produce per-cell samples by encoding control cells, adding delta_p,
    and decoding. Pair with held-out observed populations for ConfPert.

    split_type='cross_cell_line' trains on Replogle K562 (gene-intersected
    with RPE1) and tests on RPE1 ctrl input, calibrates on K562 held-out.
    """
    print(f"[scgen confpert] FUNCTION ENTRY dataset={dataset} epochs={epochs} "
          f"split={split_type}", flush=True)
    _setup_env()
    import time
    import warnings
    warnings.filterwarnings("ignore")

    import numpy as np
    import torch
    import sys
    sys.path.insert(0, "/root/src")

    # scgen 2.1.0 expects `from scvi._compat import Literal`; scvi-tools >=1.4
    # removed `_compat`. Provide a shim regardless of the installed scvi version
    # so the import doesn't break downstream of the version pin.
    import scvi
    if not hasattr(scvi, "_compat"):
        import types
        from typing import Literal as _Literal
        _shim = types.ModuleType("scvi._compat")
        _shim.Literal = _Literal
        scvi._compat = _shim
        sys.modules["scvi._compat"] = _shim
    import scgen

    print(f"[scgen confpert] CUDA={torch.cuda.is_available()}", flush=True)

    if split_type == "cross_cell_line":
        # K562 -> RPE1 transfer
        (adata, pert_col, ctrl_value, X_ctrl_calib, X_ctrl_rpe1,
         X_pert_calib_dict, X_pert_test_dict, common_perts,
         _gene_names) = _load_cross_cell_line_data()
        print(f"[scgen confpert xcl] train_adata: n={adata.shape[0]}, "
              f"d={adata.shape[1]}, common_perts={len(common_perts)}", flush=True)
    else:
        ds, adata, pert_col, is_ctrl_mask = _load_perturb_ds(
            dataset, n_top_perturbations=n_top_perturbations
        )
        print(f"[scgen confpert] adata: n={adata.shape[0]}, d={adata.shape[1]}, "
              f"n_perts={len(ds.X_pert)}, n_ctrl={int(is_ctrl_mask.sum())}", flush=True)

    # scGen needs `condition` and `cell_type` keys for its setup. Norman has no
    # meaningful cell_type partition (all K562); add a constant column.
    adata.obs["condition"] = adata.obs[pert_col].astype(str)
    adata.obs["cell_type"] = "K562"

    scgen.SCGEN.setup_anndata(adata, batch_key="condition", labels_key="cell_type")
    model = scgen.SCGEN(adata)

    print(f"[scgen confpert] training scGen for {epochs} ep ...", flush=True)
    t0 = time.time()
    model.train(
        max_epochs=epochs,
        batch_size=128,
        early_stopping=True,
        early_stopping_patience=10,
    )
    fit_t = time.time() - t0
    print(f"[scgen confpert] trained in {fit_t:.0f}s", flush=True)

    # Encode all cells to latent space directly via the underlying encoder.
    # scvi-tools 1.3 changed get_latent_representation to expect a qzv tensor
    # that scgen's encoder doesn't expose as a top-level field; bypass the
    # high-level wrapper by calling z_encoder.forward() in batches.
    model.module.eval()
    all_latent = []
    X_arr = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    X_arr = X_arr.astype(np.float32)
    batch = 256
    with torch.no_grad():
        for i in range(0, X_arr.shape[0], batch):
            chunk = torch.from_numpy(X_arr[i:i + batch]).to(model.device)
            out = model.module.z_encoder(chunk)
            qzm = out[0] if isinstance(out, tuple) else out
            all_latent.append(qzm.cpu().numpy())
    latent = np.concatenate(all_latent, axis=0)  # [n_cells, n_latent]

    labels = adata.obs[pert_col].astype(str).values
    import pandas as _pd_helper
    ctrl_value = _find_ctrl_value(_pd_helper.Series(labels))
    if ctrl_value is None:
        raise ValueError(f"no ctrl in {pert_col}")

    ctrl_mask_arr = (labels == ctrl_value)
    z_ctrl_mean = latent[ctrl_mask_arr].mean(axis=0)
    print(f"[scgen confpert] ctrl mean computed; n_ctrl_latent={int(ctrl_mask_arr.sum())}",
          flush=True)

    # Per-pert delta in latent space.
    # In within_perturbation mode, iterate over ds.X_pert (the dataset's perts).
    # In cross_cell_line mode, iterate over common_perts (the K562 ∩ RPE1 perts
    # the helper already filtered).
    if split_type == "cross_cell_line":
        pert_iter = common_perts
    else:
        pert_iter = list(ds.X_pert.keys())
    deltas = {}
    for p in pert_iter:
        mask_p = (labels == p)
        if mask_p.sum() < 5:
            continue
        deltas[p] = latent[mask_p].mean(axis=0) - z_ctrl_mean

    # Within-perturbation prediction: predict for held-out perts using K562
    # ctrl cells from the same training distribution.
    if split_type != "cross_cell_line":
        perts_to_eval = sorted(ds.X_pert.keys(),
                                key=lambda p: -ds.X_pert[p].shape[0])[:n_eval_perts]
        pred_pops, obs_pops = {}, {}
        rng = np.random.RandomState(42)
        ctrl_idx = np.where(ctrl_mask_arr)[0]
        for p in perts_to_eval:
            if p not in deltas:
                continue
            n_cells = ds.X_pert[p].shape[0]
            if len(ctrl_idx) >= n_cells:
                inp_idx = rng.choice(ctrl_idx, n_cells, replace=False)
            else:
                inp_idx = rng.choice(ctrl_idx, n_cells, replace=True)
            try:
                z_inp = latent[inp_idx]
                z_pred = z_inp + deltas[p][None, :]
                # Decode through scGen's VAE module via the canonical generative()
                # wrapper (returns dict(px=decoder(z))). Per upstream
                # github.com/theislab/scgen/blob/master/scgen/_scgenvae.py
                model.module.eval()
                with torch.no_grad():
                    z_t = torch.from_numpy(z_pred).float().to(model.device)
                    gen_out = model.module.generative(z_t)
                    px = gen_out["px"] if isinstance(gen_out, dict) else gen_out
                    X_pred = px.cpu().numpy().astype(np.float32)
                if X_pred.shape[1] != ds.X_pert[p].shape[1]:
                    print(f"  [scgen] dim mismatch {p}: pred={X_pred.shape}, "
                          f"obs={ds.X_pert[p].shape}", flush=True)
                    continue
                pred_pops[p] = X_pred
                obs_pops[p] = ds.X_pert[p]
            except Exception as e:
                print(f"  [scgen] predict failed for {p}: {type(e).__name__}: {e}",
                      flush=True)
    else:
        # Suppress the unused-variable lint for the within-pert pred_pops/obs_pops
        # in xcl mode; the xcl branch builds its own pred_calib/test_pops below.
        pred_pops, obs_pops = {}, {}

    if split_type == "cross_cell_line":
        # Predict twice per common pert: from K562 ctrl input (calib arm) and
        # from RPE1 ctrl input (test arm). Deltas were computed in latent space
        # from the K562 train split and stay the same.
        def _scgen_predict_pop(X_ctrl_input, p, n_cells):
            with torch.no_grad():
                X_arr = X_ctrl_input.astype(np.float32)
                z_inp_chunks = []
                for i in range(0, X_arr.shape[0], 256):
                    chunk = torch.from_numpy(X_arr[i:i + 256]).to(model.device)
                    out_e = model.module.z_encoder(chunk)
                    qzm = out_e[0] if isinstance(out_e, tuple) else out_e
                    z_inp_chunks.append(qzm.cpu().numpy())
                z_inp_all = np.concatenate(z_inp_chunks, axis=0)
                rng_local = np.random.RandomState(42)
                if z_inp_all.shape[0] >= n_cells:
                    sel = rng_local.choice(z_inp_all.shape[0], n_cells, replace=False)
                else:
                    sel = rng_local.choice(z_inp_all.shape[0], n_cells, replace=True)
                z_inp = z_inp_all[sel]
                z_pred = z_inp + deltas[p][None, :]
                z_t = torch.from_numpy(z_pred).float().to(model.device)
                gen_out = model.module.generative(z_t)
                px = gen_out["px"] if isinstance(gen_out, dict) else gen_out
                return px.cpu().numpy().astype(np.float32)

        pred_calib_pops, obs_calib_pops = {}, {}
        pred_test_pops, obs_test_pops = {}, {}
        for p in common_perts:
            if p not in deltas:
                continue
            try:
                n_calib = X_pert_calib_dict[p].shape[0]
                n_test = X_pert_test_dict[p].shape[0]
                pred_calib_pops[p] = _scgen_predict_pop(X_ctrl_calib, p, n_calib)
                pred_test_pops[p] = _scgen_predict_pop(X_ctrl_rpe1, p, n_test)
                obs_calib_pops[p] = X_pert_calib_dict[p]
                obs_test_pops[p] = X_pert_test_dict[p]
            except Exception as e:
                print(f"  [scgen xcl] {p}: {type(e).__name__}: {e}", flush=True)
        print(f"[scgen confpert xcl] paired {len(pred_calib_pops)} calib + "
              f"{len(pred_test_pops)} test; running conformal ...", flush=True)
        out = _run_conformal_eval_xcl(
            pred_calib_pops, obs_calib_pops, pred_test_pops, obs_test_pops,
            "replogle_k562_to_rpe1", "scgen", fit_t,
            extra_config={"epochs": epochs}
        )
        volume.commit()
        print(f"[scgen confpert xcl] -> {VOLUME_MOUNT}/confpert_scgen_xcl_"
              f"replogle_k562_to_rpe1.json", flush=True)
        return out

    print(f"[scgen confpert] paired {len(pred_pops)} perts; running conformal ...",
          flush=True)
    out = _run_conformal_eval(pred_pops, obs_pops, dataset, "scgen", fit_t,
                               extra_config={"epochs": epochs,
                                             "n_eval_perts": n_eval_perts})
    volume.commit()
    print(f"[scgen confpert] -> {VOLUME_MOUNT}/confpert_scgen_{dataset}.json",
          flush=True)
    return out


@app.local_entrypoint()
def scgen_calibrate_norman():
    """scGen + ConfPert on Norman. ~$6, ~1h on A100."""
    _scgen_calibrate_fn.remote(dataset="norman", epochs=100)


@app.local_entrypoint()
def scgen_calibrate(dataset: str = "norman", epochs: int = 100,
                     split_type: str = "within_perturbation"):
    """Parameterized scGen entrypoint.

    dataset in {norman, replogle_k562, replogle_rpe1, adamson, tahoe} for
    within_perturbation; for split_type=cross_cell_line, dataset is ignored
    and the K562->RPE1 transfer is run.
    """
    _scgen_calibrate_fn.remote(dataset=dataset, epochs=epochs,
                                 split_type=split_type)


# ---------------------------------------------------------------------------
# biolord (Piran et al. 2024 Nat Biotech; nitzanlab/biolord)
# Wrapped per lit_notes/biolord_2024.md: predict() returns (mean, var) AnnData
# pair; sample via Gaussian with calibrated noise per cell.
# ---------------------------------------------------------------------------


@app.function(
    gpu="A100",
    timeout=60 * 60 * 3,
    volumes={VOLUME_MOUNT: volume},
    image=BIOLORD_IMAGE,
)
def _biolord_calibrate_fn(dataset: str = "norman", epochs: int = 100,
                           n_eval_perts: int = 30, n_top_perturbations: int = 50,
                           split_type: str = "within_perturbation"):
    """Train biolord on a Norman-style Perturb-seq dataset, then predict
    per-cell (mean, variance) for held-out perturbations and sample Gaussians
    from those for ConfPert calibration.

    split_type='cross_cell_line' trains on Replogle K562 (gene-intersected
    with RPE1) and tests on RPE1 ctrl input, calibrates on K562 held-out.
    """
    print(f"[biolord confpert] FUNCTION ENTRY dataset={dataset} epochs={epochs} "
          f"split={split_type}", flush=True)
    _setup_env()
    import time
    import warnings
    warnings.filterwarnings("ignore")

    import numpy as np
    import torch
    import sys
    sys.path.insert(0, "/root/src")

    import biolord

    print(f"[biolord confpert] CUDA={torch.cuda.is_available()}", flush=True)

    if split_type == "cross_cell_line":
        (adata, pert_col, ctrl_value, X_ctrl_calib, X_ctrl_rpe1,
         X_pert_calib_dict, X_pert_test_dict, common_perts,
         _gene_names) = _load_cross_cell_line_data()
        print(f"[biolord confpert xcl] train_adata: n={adata.shape[0]}, "
              f"d={adata.shape[1]}, common_perts={len(common_perts)}", flush=True)
    else:
        ds, adata, pert_col, is_ctrl_mask = _load_perturb_ds(
            dataset, n_top_perturbations=n_top_perturbations
        )
        print(f"[biolord confpert] adata: n={adata.shape[0]}, d={adata.shape[1]}, "
              f"n_perts={len(ds.X_pert)}", flush=True)

    # biolord setup: categorical attribute = perturbation. Cast to str first
    # because Adamson's perturbation column has mixed types (str + NaN floats)
    # which break np.unique.sort() inside Biolord.setup_anndata at retrieval-key
    # construction (TypeError: '<' not supported between str and float).
    adata.obs[pert_col] = adata.obs[pert_col].astype(str).astype("category")

    biolord.Biolord.setup_anndata(
        adata,
        ordered_attributes_keys=None,
        categorical_attributes_keys=[pert_col],
        retrieval_attribute_key=pert_col,
    )
    model = biolord.Biolord(
        adata=adata,
        n_latent=32,
        module_params={
            "n_latent_attribute_categorical": 8,
            "gene_likelihood": "normal",
            "reconstruction_penalty": 1e2,
            "unknown_attribute_penalty": 1e1,
            "unknown_attribute_noise_param": 1e-1,
            "attribute_dropout_rate": 0.0,
            "use_batch_norm": False,
            "use_layer_norm": False,
            "seed": 42,
        },
    )

    print(f"[biolord confpert] training biolord for {epochs} ep ...", flush=True)
    t0 = time.time()
    model.train(
        max_epochs=epochs,
        batch_size=256,
        early_stopping=True,
        early_stopping_patience=10,
        plan_kwargs={"step_size_lr": 45},
    )
    fit_t = time.time() - t0
    print(f"[biolord confpert] trained in {fit_t:.0f}s", flush=True)

    labels = adata.obs[pert_col].astype(str).values
    if split_type != "cross_cell_line":
        import pandas as _pd_helper
        ctrl_value = _find_ctrl_value(_pd_helper.Series(labels))

    rng = np.random.RandomState(42)

    def _biolord_predict_pop_with_ctrl(X_ctrl_input, p, n_cells, rng_local):
        """Predict cell pop for pert p given external ctrl input cells.

        Substitutes X in a scaffold AnnData sliced from train adata so that
        biolord's setup_anndata registry is preserved while the actual cell
        expressions come from X_ctrl_input.
        """
        if X_ctrl_input.shape[0] >= n_cells:
            sel_input = rng_local.choice(X_ctrl_input.shape[0], n_cells,
                                            replace=False)
        else:
            sel_input = rng_local.choice(X_ctrl_input.shape[0], n_cells,
                                            replace=True)
        scaffold_idx = rng_local.choice(adata.shape[0], n_cells,
                                          replace=(adata.shape[0] < n_cells))
        scaffold = adata[scaffold_idx].copy()
        scaffold.X = X_ctrl_input[sel_input].astype(np.float32)
        scaffold.obs[pert_col] = p
        scaffold.obs[pert_col] = scaffold.obs[pert_col].astype(
            adata.obs[pert_col].dtype
        )
        pred_out = model.predict(adata=scaffold, nullify_attribute=[])
        if isinstance(pred_out, tuple) and len(pred_out) == 2:
            mean_ad, var_ad = pred_out
            mu = (mean_ad.X.toarray() if hasattr(mean_ad.X, "toarray")
                  else np.asarray(mean_ad.X)).astype(np.float32)
            var = (var_ad.X.toarray() if hasattr(var_ad.X, "toarray")
                   else np.asarray(var_ad.X)).astype(np.float32)
        else:
            ad_pred = pred_out
            mu = (ad_pred.X.toarray() if hasattr(ad_pred.X, "toarray")
                  else np.asarray(ad_pred.X)).astype(np.float32)
            var = np.full_like(mu, 0.05 ** 2)
        sigma = np.sqrt(np.clip(var, 1e-9, None))
        eps = rng_local.randn(*mu.shape).astype(np.float32)
        return (mu + sigma * eps).astype(np.float32)

    if split_type == "cross_cell_line":
        # Predict twice per common pert: from K562 ctrl input (calib arm) +
        # from RPE1 ctrl input (test arm).
        pred_calib_pops, obs_calib_pops = {}, {}
        pred_test_pops, obs_test_pops = {}, {}
        for p in common_perts:
            try:
                n_calib = X_pert_calib_dict[p].shape[0]
                n_test = X_pert_test_dict[p].shape[0]
                rng_p = np.random.RandomState(hash(p) % (2**31))
                pred_calib_pops[p] = _biolord_predict_pop_with_ctrl(
                    X_ctrl_calib, p, n_calib, rng_p)
                pred_test_pops[p] = _biolord_predict_pop_with_ctrl(
                    X_ctrl_rpe1, p, n_test, rng_p)
                obs_calib_pops[p] = X_pert_calib_dict[p]
                obs_test_pops[p] = X_pert_test_dict[p]
            except Exception as e:
                print(f"  [biolord xcl] {p}: {type(e).__name__}: {e}", flush=True)
        print(f"[biolord confpert xcl] paired {len(pred_calib_pops)} calib + "
              f"{len(pred_test_pops)} test; running conformal ...", flush=True)
        out = _run_conformal_eval_xcl(
            pred_calib_pops, obs_calib_pops, pred_test_pops, obs_test_pops,
            "replogle_k562_to_rpe1", "biolord", fit_t,
            extra_config={"epochs": epochs}
        )
        volume.commit()
        print(f"[biolord confpert xcl] -> {VOLUME_MOUNT}/confpert_biolord_xcl_"
              f"replogle_k562_to_rpe1.json", flush=True)
        return out

    # Within-perturbation prediction
    perts_to_eval = sorted(ds.X_pert.keys(),
                            key=lambda p: -ds.X_pert[p].shape[0])[:n_eval_perts]
    pred_pops, obs_pops = {}, {}
    ctrl_idx = np.where(is_ctrl_mask)[0]

    for p in perts_to_eval:
        n_cells = ds.X_pert[p].shape[0]
        if len(ctrl_idx) >= n_cells:
            inp_idx = rng.choice(ctrl_idx, n_cells, replace=False)
        else:
            inp_idx = rng.choice(ctrl_idx, n_cells, replace=True)
        try:
            ctrl_slice = adata[inp_idx].copy()
            ctrl_slice.obs[pert_col] = p
            ctrl_slice.obs[pert_col] = ctrl_slice.obs[pert_col].astype(
                adata.obs[pert_col].dtype
            )
            # biolord predict returns (mean_adata, var_adata)
            pred_out = model.predict(adata=ctrl_slice, nullify_attribute=[])
            if isinstance(pred_out, tuple) and len(pred_out) == 2:
                mean_ad, var_ad = pred_out
                mu = (mean_ad.X.toarray() if hasattr(mean_ad.X, "toarray")
                      else np.asarray(mean_ad.X)).astype(np.float32)
                var = (var_ad.X.toarray() if hasattr(var_ad.X, "toarray")
                       else np.asarray(var_ad.X)).astype(np.float32)
            else:
                # Fallback: API returned only the mean
                ad = pred_out
                mu = (ad.X.toarray() if hasattr(ad.X, "toarray")
                      else np.asarray(ad.X)).astype(np.float32)
                var = np.full_like(mu, 0.05 ** 2)
            sigma = np.sqrt(np.clip(var, 1e-9, None))
            eps = rng.randn(*mu.shape).astype(np.float32)
            X_pred = (mu + sigma * eps).astype(np.float32)
            if X_pred.shape[1] != ds.X_pert[p].shape[1]:
                print(f"  [biolord] dim mismatch {p}: pred={X_pred.shape}, "
                      f"obs={ds.X_pert[p].shape}", flush=True)
                continue
            pred_pops[p] = X_pred
            obs_pops[p] = ds.X_pert[p]
        except Exception as e:
            print(f"  [biolord] predict failed for {p}: {type(e).__name__}: {e}",
                  flush=True)

    print(f"[biolord confpert] paired {len(pred_pops)} perts; running conformal ...",
          flush=True)
    out = _run_conformal_eval(pred_pops, obs_pops, dataset, "biolord", fit_t,
                               extra_config={"epochs": epochs,
                                             "n_eval_perts": n_eval_perts})
    volume.commit()
    print(f"[biolord confpert] -> {VOLUME_MOUNT}/confpert_biolord_{dataset}.json",
          flush=True)
    return out


@app.local_entrypoint()
def biolord_calibrate_norman():
    """biolord + ConfPert on Norman. ~$36, ~2h on A100."""
    _biolord_calibrate_fn.remote(dataset="norman", epochs=100)


@app.local_entrypoint()
def biolord_calibrate(dataset: str = "norman", epochs: int = 100,
                       split_type: str = "within_perturbation"):
    """Parameterized biolord entrypoint.

    dataset in {norman, replogle_k562, replogle_rpe1, adamson, tahoe} for
    within_perturbation; for split_type=cross_cell_line, dataset is ignored.
    """
    _biolord_calibrate_fn.remote(dataset=dataset, epochs=epochs,
                                   split_type=split_type)


# ---------------------------------------------------------------------------
# sVAE+ / SAMS-VAE (Bereket & Karaletsos 2023 NeurIPS; insitro/sams-vae)
# Wrapped per lit_notes/svaeplus_2023.md: per-cell posterior samples via
# SAMSVAEPredictor.sample_observations(). Build SAMSVAEModel + correlated
# normal guide directly to bypass the YAML-driven Hydra training.
# ---------------------------------------------------------------------------


@app.function(
    gpu="A100",
    timeout=60 * 60 * 3,
    volumes={VOLUME_MOUNT: volume},
    image=SVAEPLUS_IMAGE,
)
def _svaeplus_calibrate_fn(dataset: str = "norman", epochs: int = 30,
                            n_eval_perts: int = 30, n_top_perturbations: int = 50,
                            n_latent: int = 32, batch_size: int = 128):
    """Train SAMS-VAE on a Norman-style dataset, then sample per-cell
    counterfactuals via SAMSVAEPredictor.sample_observations() and feed
    (pred_pop, obs_pop) pairs into ConfPert's PerturbationConformal head.
    """
    print(f"[svae+ confpert] FUNCTION ENTRY dataset={dataset} epochs={epochs}",
          flush=True)
    _setup_env()
    import time
    import warnings
    warnings.filterwarnings("ignore")

    import numpy as np
    import torch
    import sys
    sys.path.insert(0, "/root/src")

    from sams_vae.models.sams_vae.model import SAMSVAEModel
    from sams_vae.models.sams_vae.guides.correlated_normal_guide import (
        SAMSVAECorrelatedNormalGuide,
    )
    from sams_vae.models.sams_vae.loss_module import SAMSVAE_ELBOLossModule
    from sams_vae.models.sams_vae.predictor import SAMSVAEPredictor
    from sams_vae.data.utils.perturbation_dataset import (
        SCRNASeqTensorPerturbationDataset,
    )
    from sams_vae.data.utils.perturbation_datamodule import (
        ObservationNormalizationStatistics,
    )

    print(f"[svae+ confpert] CUDA={torch.cuda.is_available()}", flush=True)

    ds, adata, pert_col, is_ctrl_mask = _load_perturb_ds(
        dataset, n_top_perturbations=n_top_perturbations
    )
    print(f"[svae+ confpert] adata: n={adata.shape[0]}, d={adata.shape[1]}, "
          f"n_perts={len(ds.X_pert)}, n_ctrl={int(is_ctrl_mask.sum())}", flush=True)

    # Build one-hot perturbation encoding for SAMS-VAE's D tensor
    labels = adata.obs[pert_col].astype(str).values
    unique_perts = sorted(set(labels))
    pert_to_idx = {p: i for i, p in enumerate(unique_perts)}
    n_treatments = len(unique_perts)
    n_phenos = adata.shape[1]

    X_arr = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    X_arr = X_arr.astype(np.float32)
    D_arr = np.zeros((adata.shape[0], n_treatments), dtype=np.float32)
    for i, p in enumerate(labels):
        D_arr[i, pert_to_idx[p]] = 1.0
    X_t = torch.from_numpy(X_arr)
    D_t = torch.from_numpy(D_arr)

    # 80/20 train/val split (held-out perts come from K1 calibration loop later)
    rng = np.random.RandomState(0)
    idx = np.arange(adata.shape[0])
    rng.shuffle(idx)
    cut = int(0.8 * adata.shape[0])
    tr_idx, vl_idx = idx[:cut], idx[cut:]

    train_dataset = SCRNASeqTensorPerturbationDataset(X=X_t[tr_idx], D=D_t[tr_idx])

    # Normalization stats per their convention
    x_train = X_t[tr_idx]
    x_norm_stats = ObservationNormalizationStatistics(
        x_mean=x_train.mean(0),
        x_std=x_train.std(0).clamp(min=1e-6),
        log_x_mean=torch.log1p(x_train).mean(0),
        log_x_std=torch.log1p(x_train).std(0).clamp(min=1e-6),
    )

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = SAMSVAEModel(
        n_latent=n_latent,
        n_treatments=n_treatments,
        n_phenos=n_phenos,
        mask_prior_prob=0.01,
        embedding_prior_scale=1.0,
        likelihood_key="normal",
        decoder_n_layers=2,
        decoder_n_hidden=256,
    ).to(device)
    guide = SAMSVAECorrelatedNormalGuide(
        n_latent=n_latent,
        n_treatments=n_treatments,
        n_phenos=n_phenos,
        basal_encoder_n_layers=2,
        basal_encoder_n_hidden=256,
        # Our X is already log1p-normalized via _load_perturb_ds; use plain
        # standardize (z-score with x_mean / x_std), not log_standardize
        # (which re-applies log1p and produces NaN on already-logged inputs).
        basal_encoder_input_normalization="standardize",
        embedding_encoder_n_layers=2,
        embedding_encoder_n_hidden=256,
        x_normalization_stats=x_norm_stats,
        gs_temperature=1.0,
    ).to(device)

    loss_module = SAMSVAE_ELBOLossModule(model=model, guide=guide)
    # Lower LR + grad clipping for stability (1e-3 + no-clip diverged to NaN
    # in v5; the SAMS-VAE paper uses lr=5e-4 in their reference YAML).
    optim = torch.optim.Adam(
        list(model.parameters()) + list(guide.parameters()), lr=5e-4
    )

    train_loader = torch.utils.data.DataLoader(
        train_dataset, batch_size=batch_size, shuffle=True
    )
    n_train = len(train_dataset)
    # D_obs_counts: per-perturbation total observation count, used by the
    # ELBO loss to reweight perturbation-plated variables. shape (n_treatments,)
    D_obs_counts = (D_t[tr_idx] != 0).sum(0).float().to(device)

    print(f"[svae+ confpert] training SAMS-VAE for {epochs} epochs, "
          f"n_train={n_train}, batches/epoch={len(train_loader)} ...", flush=True)
    t0 = time.time()
    for epoch in range(epochs):
        epoch_loss = 0.0
        n_batches = 0
        for batch in train_loader:
            X_b = batch["X"].to(device)
            D_b = batch["D"].to(device)
            optim.zero_grad()
            # PerturbationPlatedELBOLossModule.loss(X, D, D_obs_counts, ...)
            # returns (loss, metrics_dict). Module.forward returns the lower-level
            # (guide_dists, model_dists, samples) tuple, NOT a backprop-able loss.
            loss, _metrics = loss_module.loss(
                X=X_b, D=D_b, D_obs_counts=D_obs_counts, n_particles=1
            )
            loss.backward()
            # Grad clip @ 5.0 prevents the divergence-to-NaN observed in v5
            # (mask logits exploded once the ELBO went strongly negative).
            torch.nn.utils.clip_grad_norm_(
                list(model.parameters()) + list(guide.parameters()), max_norm=5.0
            )
            optim.step()
            epoch_loss += float(loss.detach())
            n_batches += 1
        if (epoch + 1) % 5 == 0 or epoch == 0:
            print(f"[svae+ confpert] epoch {epoch + 1}/{epochs} "
                  f"loss={epoch_loss / max(n_batches, 1):.3f}", flush=True)
    fit_t = time.time() - t0
    print(f"[svae+ confpert] trained in {fit_t:.0f}s", flush=True)

    # Build per-pert (pred, obs) pairs via SAMSVAEPredictor.sample_observations
    predictor = SAMSVAEPredictor(model=model, guide=guide).to(device)

    pred_pops, obs_pops = {}, {}
    perts_to_eval = sorted(ds.X_pert.keys(),
                           key=lambda p: -ds.X_pert[p].shape[0])[:n_eval_perts]
    for p in perts_to_eval:
        if p not in pert_to_idx:
            continue
        n_cells = ds.X_pert[p].shape[0]
        # Build dosage tensor: one-hot for the target pert, repeated n_cells times
        D_pert = torch.zeros((1, n_treatments), dtype=torch.float32, device=device)
        D_pert[0, pert_to_idx[p]] = 1.0
        try:
            sample_adata = predictor.sample_observations(
                dosages=D_pert, perturbation_names=[p], n_particles=n_cells
            )
            X_pred = np.asarray(sample_adata.X, dtype=np.float32)
            if X_pred.shape[0] != n_cells:
                # repeat to match observed count if undersampled
                k = (n_cells + X_pred.shape[0] - 1) // X_pred.shape[0]
                X_pred = np.tile(X_pred, (k, 1))[:n_cells]
            if X_pred.shape[1] != ds.X_pert[p].shape[1]:
                print(f"  [svae+] dim mismatch {p}: pred={X_pred.shape}, "
                      f"obs={ds.X_pert[p].shape}", flush=True)
                continue
            pred_pops[p] = X_pred
            obs_pops[p] = ds.X_pert[p]
        except Exception as e:
            print(f"  [svae+] predict failed for {p}: {type(e).__name__}: {e}",
                  flush=True)

    print(f"[svae+ confpert] paired {len(pred_pops)} perts; running conformal ...",
          flush=True)
    out = _run_conformal_eval(pred_pops, obs_pops, dataset, "svaeplus", fit_t,
                               extra_config={"epochs": epochs,
                                             "n_eval_perts": n_eval_perts,
                                             "n_latent": n_latent})
    volume.commit()
    print(f"[svae+ confpert] -> {VOLUME_MOUNT}/confpert_svaeplus_{dataset}.json",
          flush=True)
    return out


@app.local_entrypoint()
def svaeplus_calibrate(dataset: str = "norman", epochs: int = 30):
    """Parameterized sVAE+ / SAMS-VAE entrypoint."""
    _svaeplus_calibrate_fn.remote(dataset=dataset, epochs=epochs)


# ---------------------------------------------------------------------------
# Tahoe-100M subset builder (Vevo / Arc VCA 2025; tahoebio/Tahoe-100M).
#
# Streams the 429 GB HF parquet shards, filters to drugs that overlap PRISM
# 24Q2 + N representative cell lines, materializes one h5ad subset
# (target ~500 K cells per compute_estimate.md) on the causeflow-artifacts
# volume. CPU-only: filtering is I/O-bound, no GPU needed.
#
# Prereq: PRISM 24Q2 extended_compound_list.csv must already be on the
# volume at /artifacts/data/k3/prism_24q2_extended_compound_list.csv.
# `scripts/download_k3_data.py` already places it there during K3 setup.
# ---------------------------------------------------------------------------

TAHOE_BUILDER_IMAGE = _add_local(
    _BASE_DEPS.pip_install(
        "datasets>=2.20",
        "huggingface_hub>=0.20",
        "pyarrow>=14",
        "fsspec>=2024.2",
    )
)


@app.function(
    cpu=4.0,
    memory=24 * 1024,  # 24 GB; AnnData materialization is the hot spot
    timeout=60 * 60 * 8,  # 8 h: streaming-filter is I/O-bound
    volumes={VOLUME_MOUNT: volume},
    image=TAHOE_BUILDER_IMAGE,
    secrets=[modal.Secret.from_name("huggingface")],
)
def _tahoe_subset_build_fn(
    n_drugs: int = 50,
    n_cell_lines: int = 10,
    cells_per_condition: int = 400,
    dose_um: float = 5.0,
    output_basename: str = "tahoe_subset",
    prism_compounds_remote_csv: str = "/artifacts/data/k3/prism_24q2_extended_compound_list.csv",
    max_shards: int | None = None,
    max_stagnant_shards: int = 200,
) -> dict:
    """Build a Tahoe-100M subset h5ad on the Modal volume.

    Selection rule:
      1. Compute the case-insensitive intersection of Tahoe ``drug_metadata.drug``
         names with PRISM 24Q2 ``Drug.Name``. (PRISM Drug.Name is curated;
         Tahoe ``drug`` is the trade / generic name.)
      2. Pick top-``n_cell_lines`` Tahoe cell lines by total cell count among
         those drugs at the chosen dose. (Cellosaurus IDs e.g. ``CVCL_0023``.)
      3. Stream the main expression_data, keep rows where
         ``drug in matched_drugs``, ``cell_line_id in chosen_lines``, and
         (for non-DMSO rows) the dose annotation matches ``dose_um``.
      4. Subsample to at most ``cells_per_condition`` cells per
         ``(drug, cell_line_id)`` to balance the subset.
      5. Always retain DMSO_TF rows (Tahoe vehicle control) for the chosen
         cell lines, since they are the K1 / K3 control population.
      6. Materialize an AnnData with raw counts in ``X``, drug + cell line in
         ``obs``, gene symbols (mapped via ``gene_metadata`` from token IDs)
         in ``var_names``, save as h5ad on the volume.

    Returns a manifest dict. Output paths:
      /artifacts/data/{output_basename}.h5ad
      /artifacts/data/{output_basename}_manifest.json
    """
    print(f"[tahoe-subset] FUNCTION ENTRY n_drugs={n_drugs} "
          f"n_cell_lines={n_cell_lines} dose_um={dose_um}", flush=True)
    _setup_env()
    import json
    import os
    import time
    import warnings
    warnings.filterwarnings("ignore")

    import numpy as np
    import pandas as pd
    import anndata as ad
    from datasets import load_dataset
    from scipy import sparse

    # Re-export the HF token from the modal secret so `datasets` finds it.
    hf_token = (os.environ.get("HF_TOKEN")
                or os.environ.get("HUGGING_FACE_HUB_TOKEN"))
    if hf_token:
        os.environ["HF_TOKEN"] = hf_token
        os.environ["HUGGING_FACE_HUB_TOKEN"] = hf_token

    out_dir = f"{VOLUME_MOUNT}/data"
    os.makedirs(out_dir, exist_ok=True)
    out_h5ad = f"{out_dir}/{output_basename}.h5ad"
    out_manifest = f"{out_dir}/{output_basename}_manifest.json"

    # ------------------------------------------------------------------
    # Step 1: load PRISM compound list (already on volume), Tahoe drug+cell
    # metadata configs (small tables, fully materialize).
    # ------------------------------------------------------------------
    def _canon_drug(s: str) -> str:
        """Canonicalise drug names for cross-DB matching: uppercase, strip
        whitespace, hyphens, parentheses, plus/dot. PRISM uses ``AMG-397``;
        Tahoe may use ``AMG397`` or ``AMG-397`` etc. Strip everything
        non-alphanumeric so the comparison is stable.
        """
        import re
        return re.sub(r"[^A-Z0-9]+", "", str(s).upper())

    if not os.path.exists(prism_compounds_remote_csv):
        raise FileNotFoundError(
            f"PRISM compound CSV not found at {prism_compounds_remote_csv}. "
            "Run scripts/download_k3_data.py first to populate the volume."
        )
    prism = pd.read_csv(prism_compounds_remote_csv)
    prism_names_canon = set(prism["Drug.Name"].astype(str).map(_canon_drug))
    prism_names_canon.discard("")
    prism_names_canon.discard("DMSO")
    print(f"[tahoe-subset] {len(prism_names_canon)} unique PRISM drug names "
          f"(canonicalised)", flush=True)

    print("[tahoe-subset] loading Tahoe drug_metadata config ...", flush=True)
    drug_meta = load_dataset(
        "tahoebio/Tahoe-100M", "drug_metadata", split="train", token=hf_token,
    ).to_pandas()
    print(f"[tahoe-subset]   drug_metadata rows={len(drug_meta)} "
          f"cols={list(drug_meta.columns)}", flush=True)

    print("[tahoe-subset] loading Tahoe cell_line_metadata config ...",
          flush=True)
    cell_meta = load_dataset(
        "tahoebio/Tahoe-100M", "cell_line_metadata", split="train", token=hf_token,
    ).to_pandas()
    print(f"[tahoe-subset]   cell_line_metadata rows={len(cell_meta)} "
          f"cols={list(cell_meta.columns)}", flush=True)

    print("[tahoe-subset] loading Tahoe gene_metadata config ...", flush=True)
    gene_meta = load_dataset(
        "tahoebio/Tahoe-100M", "gene_metadata", split="train", token=hf_token,
    ).to_pandas()
    print(f"[tahoe-subset]   gene_metadata rows={len(gene_meta)} "
          f"cols={list(gene_meta.columns)}", flush=True)

    # ------------------------------------------------------------------
    # Step 2: drug overlap (canonicalised case- and punctuation-insensitive
    # match). Strip dashes / whitespace / parens on both sides to catch the
    # common AMG-397 vs AMG397 vs AMG 397 mismatch the lit-side notes.
    # ------------------------------------------------------------------
    drug_col = ("drug" if "drug" in drug_meta.columns
                else ("name" if "name" in drug_meta.columns
                      else drug_meta.columns[0]))
    drug_meta["__canon"] = drug_meta[drug_col].astype(str).map(_canon_drug)
    overlap_mask = drug_meta["__canon"].isin(prism_names_canon)
    overlap_drugs_full = drug_meta.loc[overlap_mask, drug_col].astype(str).tolist()
    print(f"[tahoe-subset] {len(overlap_drugs_full)} Tahoe drugs match PRISM "
          f"(canonicalised match)", flush=True)
    if len(overlap_drugs_full) < 5:
        sample_tahoe = drug_meta[drug_col].head(20).tolist()
        sample_prism = list(prism_names_canon)[:20]
        raise RuntimeError(
            f"Only {len(overlap_drugs_full)} drug matches after canonicalisation. "
            f"Sample Tahoe drugs: {sample_tahoe!r}. Sample PRISM canon: "
            f"{sample_prism!r}."
        )
    overlap_drugs = set(overlap_drugs_full)

    # ------------------------------------------------------------------
    # Step 3: pick cell lines from cell_line_metadata. Prefer lines with a
    # DepMap ID (PRISM-mappable) and adequate Tahoe coverage. We do NOT
    # scan expression_data to count cells per line — Tahoe is balanced by
    # design (every plate carries the same 50-line pool, per the lit note),
    # so cell counts per (drug, line) are roughly uniform for any drug
    # actually run on a given plate.
    # ------------------------------------------------------------------
    cl_id_col = next((c for c in ("Cell_ID_Cellosaur", "Cell_ID_Cellosaurus",
                                    "cell_line_id", "Cellosaurus_ID")
                       if c in cell_meta.columns), None)
    depmap_col = next((c for c in ("Cell_ID_DepMap", "DepMap_ID", "depmap_id")
                        if c in cell_meta.columns), None)
    if cl_id_col is None:
        raise RuntimeError(
            f"cell_line_metadata has no Cellosaurus-ID column; cols="
            f"{list(cell_meta.columns)}"
        )
    cell_meta_pool = cell_meta.copy()
    if depmap_col is not None:
        # Strict: only keep DepMap-mapped lines so we can transfer to PRISM.
        cell_meta_pool = cell_meta_pool[
            cell_meta_pool[depmap_col].notna() &
            (cell_meta_pool[depmap_col].astype(str).str.len() > 0)
        ]
    # cell_line_metadata may carry multiple rows per Cellosaurus ID (one per
    # sample/replicate); deduplicate before the head() pick.
    cell_meta_pool = cell_meta_pool.drop_duplicates(subset=[cl_id_col])
    chosen_lines = (cell_meta_pool[cl_id_col].astype(str)
                                              .head(n_cell_lines)
                                              .tolist())
    chosen_lines_set = set(chosen_lines)
    print(f"[tahoe-subset] chosen_lines (n={len(chosen_lines)}, pre-stream): "
          f"{chosen_lines}", flush=True)
    if not chosen_lines:
        raise RuntimeError("zero cell lines selected; check metadata schema")

    # ------------------------------------------------------------------
    # Step 4: SINGLE streaming pass over expression_data with all filters
    # applied. PIN the config to "expression_data" — without the config arg,
    # `datasets.load_dataset` may pick the wrong default config (per-cell
    # obs_metadata has no expressions and is a silent failure mode).
    # ------------------------------------------------------------------
    # Quotas: one bucket per (drug, line); DMSO bucket is 4x per line.
    quota: dict[tuple[str, str], int] = {
        (d, cl): cells_per_condition
        for d in overlap_drugs for cl in chosen_lines_set
    }
    for cl in chosen_lines_set:
        quota[("__DMSO__", cl)] = 4 * cells_per_condition
    print(f"[tahoe-subset] {len(quota)} buckets, target ~"
          f"{sum(quota.values()):,} cells", flush=True)

    # Use pyarrow.dataset with HuggingFace's HfFileSystem; this enables
    # parquet predicate pushdown (skip whole row groups whose `drug` column
    # statistics don't intersect overlap_drugs). HF datasets streaming was
    # observed to stall on the 429 GB Tahoe expression_data — the pyarrow
    # path is dramatically faster because it filters at row-group granularity
    # and only downloads matching shards' relevant columns.
    print("[tahoe-subset] using pyarrow.dataset on HF Hub parquet shards ...",
          flush=True)
    import pyarrow as pa
    import pyarrow.dataset as pads
    import pyarrow.compute as pc
    from huggingface_hub import HfFileSystem
    fs = HfFileSystem(token=hf_token)
    # Tahoe-100M ships 3388 parquet shards (~80 MB each) on the main branch
    # at datasets/tahoebio/Tahoe-100M/data/train-NNNNN-of-03388.parquet.
    parquet_root = "datasets/tahoebio/Tahoe-100M/data"
    try:
        all_files = fs.ls(parquet_root, detail=False)
        sample_files = all_files[:3]
        n_shards = len(all_files)
        print(f"[tahoe-subset]   {n_shards} parquet shards under "
              f"{parquet_root}; sample: {sample_files}", flush=True)
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(
            f"HF parquet root {parquet_root} not accessible: {e}. "
            "Check token + dataset license."
        )

    # Iterate parquet shards individually with per-shard predicate pushdown.
    # This lets us print progress, stop early when all quotas are met, and
    # avoid materialising 20+ GB of (drug, cell_line)-matching rows in
    # memory. We process the shards in random order so a random subset of
    # the dataset gives proportional coverage across plates / cell lines /
    # drugs.
    cols_wanted = ["genes", "expressions", "drug", "cell_line_id",
                    "moa-fine", "canonical_smiles", "pubchem_cid", "plate",
                    "sample"]
    # Probe a single shard's schema to filter cols_wanted to only those that
    # actually exist (Tahoe-100M's main-branch parquet may not carry every
    # field we name).
    _probe = pads.dataset(all_files[0], filesystem=fs, format="parquet")
    cols_needed = [c for c in cols_wanted if c in _probe.schema.names]
    drug_filter = pc.is_in(pc.field("drug"),
                            value_set=pa.array(list(overlap_drugs)))
    dmso_filter = pc.starts_with(pc.field("drug"), "DMSO")
    line_filter = pc.is_in(pc.field("cell_line_id"),
                            value_set=pa.array(list(chosen_lines_set)))
    full_filter = pc.and_kleene(line_filter,
                                  pc.or_kleene(drug_filter, dmso_filter))
    print(f"[tahoe-subset]   per-shard scan with pushdown on "
          f"{len(cols_needed)} cols ...", flush=True)

    rows_genes: list = []
    rows_expr: list = []
    rows_obs: list = []
    rng_shuffle = np.random.RandomState(0)
    shard_order = list(all_files)
    rng_shuffle.shuffle(shard_order)
    if max_shards is not None:
        shard_order = shard_order[:max_shards]
        print(f"[tahoe-subset]   capped to first {len(shard_order)} shards "
              f"(max_shards={max_shards})", flush=True)
    quota_remaining = dict(quota)
    seen = 0
    kept = 0
    last_kept_change_shard = 0
    last_kept_value = 0
    t0 = time.time()
    for shard_idx, shard_path in enumerate(shard_order):
        try:
            shard_ds = pads.dataset(shard_path, filesystem=fs,
                                     format="parquet")
            shard_table = shard_ds.to_table(columns=cols_needed,
                                              filter=full_filter)
        except Exception as e:  # noqa: BLE001
            print(f"[tahoe-subset]   shard {shard_idx} read failed: {e}",
                  flush=True)
            continue
        seen += shard_table.num_rows  # rows AFTER pushdown filter
        if shard_table.num_rows == 0:
            if shard_idx % 50 == 0:
                print(f"[tahoe-subset]   shard {shard_idx}/{len(shard_order)} "
                      f"empty (elapsed {time.time()-t0:.0f}s, kept={kept:,})",
                      flush=True)
            continue
        drugs_arr = shard_table["drug"].to_numpy(zero_copy_only=False)
        lines_arr = shard_table["cell_line_id"].to_numpy(zero_copy_only=False)
        genes_col = shard_table["genes"]
        expr_col = shard_table["expressions"]
        moa_arr = (shard_table["moa-fine"].to_numpy(zero_copy_only=False)
                     if "moa-fine" in cols_needed else None)
        smi_arr = (shard_table["canonical_smiles"].to_numpy(zero_copy_only=False)
                     if "canonical_smiles" in cols_needed else None)
        pcid_arr = (shard_table["pubchem_cid"].to_numpy(zero_copy_only=False)
                      if "pubchem_cid" in cols_needed else None)
        plate_arr = (shard_table["plate"].to_numpy(zero_copy_only=False)
                       if "plate" in cols_needed else None)
        sample_arr = (shard_table["sample"].to_numpy(zero_copy_only=False)
                        if "sample" in cols_needed else None)
        for i in range(shard_table.num_rows):
            d = str(drugs_arr[i])
            cl = str(lines_arr[i])
            is_dmso = d.upper().startswith("DMSO")
            key = ("__DMSO__", cl) if is_dmso else (d, cl)
            if quota_remaining.get(key, 0) <= 0:
                continue
            quota_remaining[key] -= 1
            rows_genes.append(np.asarray(genes_col[i].as_py(), dtype=np.int64))
            rows_expr.append(np.asarray(expr_col[i].as_py(), dtype=np.float32))
            rows_obs.append({
                "drug": "DMSO" if is_dmso else d,
                "cell_line_id": cl,
                "moa_fine": str(moa_arr[i]) if moa_arr is not None else "",
                "canonical_smiles": str(smi_arr[i]) if smi_arr is not None else "",
                "pubchem_cid": str(pcid_arr[i]) if pcid_arr is not None else "",
                "plate": str(plate_arr[i]) if plate_arr is not None else "",
                "sample": str(sample_arr[i]) if sample_arr is not None else "",
            })
            kept += 1
        if shard_idx % 25 == 0 or shard_idx < 5:
            print(f"[tahoe-subset]   shard {shard_idx}/{len(shard_order)} "
                  f"matched {shard_table.num_rows:,} rows; "
                  f"total kept={kept:,} elapsed={time.time()-t0:.0f}s",
                  flush=True)
        n_unmet = sum(1 for q in quota_remaining.values() if q > 0)
        if n_unmet == 0:
            print(f"[tahoe-subset]   all quotas met after shard {shard_idx}",
                  flush=True)
            break
        # Stagnation early-exit: if we've gone N consecutive shards without
        # the kept count changing, the remaining buckets won't fill (likely
        # rare drug-line combos with little Tahoe coverage). Stop.
        if kept != last_kept_value:
            last_kept_change_shard = shard_idx
            last_kept_value = kept
        elif shard_idx - last_kept_change_shard >= max_stagnant_shards:
            print(f"[tahoe-subset]   stagnant for {max_stagnant_shards} "
                  f"shards at kept={kept:,}; stopping early. "
                  f"{n_unmet} buckets still unfilled.", flush=True)
            break
    print(f"[tahoe-subset] scan done: shards_processed<={shard_idx + 1} "
          f"matched_rows={seen:,} kept={kept:,} "
          f"({time.time()-t0:.0f}s)", flush=True)

    # Post-stream: drop any drugs we collected fewer than min_cells_per_drug
    # cells for. n_drugs is honoured by ranking the surviving drugs by total
    # cells.
    pass1_counts: dict[tuple[str, str], int] = {}
    pass1_seen = seen      # alias for backwards-compat manifest fields
    pass1_kept = kept
    pass2_seen = seen
    pass2_kept = kept
    for o in rows_obs:
        d = o["drug"]
        cl = o["cell_line_id"]
        if d == "DMSO":
            continue
        pass1_counts[(d, cl)] = pass1_counts.get((d, cl), 0) + 1
    drug_total: dict[str, int] = {}
    for (d, cl), c in pass1_counts.items():
        drug_total[d] = drug_total.get(d, 0) + c
    chosen_drugs = sorted(drug_total, key=drug_total.get, reverse=True)[:n_drugs]
    chosen_drugs_set = set(chosen_drugs)
    keep_idx = [i for i, o in enumerate(rows_obs)
                 if o["drug"] == "DMSO" or o["drug"] in chosen_drugs_set]
    rows_genes = [rows_genes[i] for i in keep_idx]
    rows_expr = [rows_expr[i] for i in keep_idx]
    rows_obs = [rows_obs[i] for i in keep_idx]
    print(f"[tahoe-subset] post-rank: kept {len(rows_obs):,} cells across "
          f"{len(chosen_drugs)} drugs + DMSO", flush=True)

    if pass2_kept == 0:
        raise RuntimeError("pass 2 collected zero cells; selection rule broken")

    # ------------------------------------------------------------------
    # Step 5: assemble AnnData. genes are token IDs into gene_meta; build
    # the full sparse cells x genes matrix. Use the full token-ID universe
    # as columns (sparse zeros elsewhere), then drop empty columns at the
    # end.
    # ------------------------------------------------------------------
    print("[tahoe-subset] assembling AnnData ...", flush=True)
    n_cells = len(rows_expr)
    indptr = [0]
    indices_list = []
    data_list = []
    n_genes_full = int(gene_meta["token_id"].max()) + 1
    for g, e in zip(rows_genes, rows_expr):
        indices_list.append(g)
        data_list.append(e)
        indptr.append(indptr[-1] + len(g))
    indices_arr = np.concatenate(indices_list, axis=0).astype(np.int32)
    data_arr = np.concatenate(data_list, axis=0).astype(np.float32)
    indptr_arr = np.asarray(indptr, dtype=np.int64)
    # Defensive: any token ID outside gene_meta's universe means we'd silently
    # write to the wrong column. Fail fast with the bad ID for debugging.
    max_token = int(indices_arr.max()) if indices_arr.size else -1
    if max_token >= n_genes_full:
        raise RuntimeError(
            f"expression token ID {max_token} >= gene_meta universe "
            f"{n_genes_full}; gene_meta is stale relative to expression_data."
        )
    X_csr = sparse.csr_matrix((data_arr, indices_arr, indptr_arr),
                               shape=(n_cells, n_genes_full))
    # Drop all-zero columns to shrink (keeps gene mapping intact via var)
    col_nnz = np.asarray(X_csr.getnnz(axis=0)).ravel()
    keep_cols = np.where(col_nnz > 0)[0]
    X_csr = X_csr[:, keep_cols]
    print(f"[tahoe-subset]   X shape after dropping zero cols: {X_csr.shape}",
          flush=True)

    # Map token_id -> gene_symbol / ensembl_id. ``reindex`` fills missing
    # token IDs with NaN instead of raising, which is what we want here:
    # we already failed-fast on out-of-universe IDs above, so any remaining
    # gaps are gene_meta rows missing a symbol — fall back to the token id.
    gene_meta_idx = gene_meta.set_index("token_id")
    var = pd.DataFrame(index=keep_cols.astype(str))
    var["token_id"] = keep_cols
    if "gene_symbol" in gene_meta_idx.columns:
        var["gene_symbol"] = (gene_meta_idx.reindex(keep_cols)["gene_symbol"]
                                            .astype(str).values)
    else:
        var["gene_symbol"] = var.index
    if "ensembl_id" in gene_meta_idx.columns:
        var["ensembl_id"] = (gene_meta_idx.reindex(keep_cols)["ensembl_id"]
                                            .astype(str).values)
    else:
        var["ensembl_id"] = ""
    # Use gene_symbol as var_names (with token_id fallback for missing)
    sym = var["gene_symbol"].astype(str)
    sym = sym.where(sym.str.len() > 0, var["token_id"].astype(str))
    # Drop the existing gene_symbol column so reset_index can promote the
    # gene_symbol-named index to a column without conflict.
    var = var.drop(columns=["gene_symbol"])
    var.index = pd.Index(sym, name="gene_symbol")
    # Dedupe by first occurrence of gene_symbol (preserving original order
    # so X_csr column alignment can be reproduced).
    keep_first = ~var.index.duplicated(keep="first")
    n_before_dedup = len(var)
    var = var[keep_first].copy()
    if len(var) != X_csr.shape[1]:
        keep_idx = np.where(keep_first)[0]
        X_csr = X_csr[:, keep_idx]
        print(f"[tahoe-subset]   deduped X {n_before_dedup}->{X_csr.shape[1]} cols",
              flush=True)

    obs = pd.DataFrame(rows_obs)
    obs.index = obs.index.astype(str)
    adata = ad.AnnData(X=X_csr, obs=obs, var=var)
    adata.uns["tahoe_subset_version"] = "v1"
    adata.uns["chosen_drugs"] = list(chosen_drugs)
    adata.uns["chosen_cell_lines"] = list(chosen_lines)
    adata.uns["dose_um"] = float(dose_um)
    adata.uns["cells_per_condition"] = int(cells_per_condition)

    print(f"[tahoe-subset] writing {out_h5ad} ({adata.shape}) ...", flush=True)
    adata.write_h5ad(out_h5ad, compression="gzip")

    manifest = {
        "out_h5ad": out_h5ad,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "chosen_drugs": list(chosen_drugs),
        "chosen_cell_lines": list(chosen_lines),
        "n_pass1_seen": int(pass1_seen),
        "n_pass1_kept": int(pass1_kept),
        "n_pass2_seen": int(pass2_seen),
        "n_pass2_kept": int(pass2_kept),
        "dose_um": float(dose_um),
        "cells_per_condition_target": int(cells_per_condition),
        "n_prism_overlap_drugs": int(len(overlap_drugs_full)),
    }
    with open(out_manifest, "w") as f:
        json.dump(manifest, f, indent=2)
    volume.commit()
    print(f"[tahoe-subset] manifest -> {out_manifest}", flush=True)
    return manifest


@app.local_entrypoint()
def tahoe_subset_build(n_drugs: int = 50, n_cell_lines: int = 10,
                        cells_per_condition: int = 400, dose_um: float = 5.0,
                        max_shards: int = 0,
                        max_stagnant_shards: int = 200):
    """Build the Tahoe-100M PRISM-overlap subset on the Modal volume.

    One-time job. Output: /artifacts/data/tahoe_subset.h5ad. After it lands,
    pull locally with `modal volume get causeflow-artifacts /data/tahoe_subset.h5ad
    /Volumes/STORAGE/FLM4S_data/tahoe_subset.h5ad` and symlink into
    `confpert/data/`.

    Tuning:
      - cells_per_condition controls bucket size (larger → more data but
        longer scan if rare buckets never fill).
      - max_shards forces a hard cap on the scan (0 = no cap, scan all 3388).
      - max_stagnant_shards triggers early-stop if the kept count hasn't
        changed in this many consecutive shards (default 200 ≈ 5 min idle).
    """
    ms = max_shards if max_shards > 0 else None
    res = _tahoe_subset_build_fn.remote(
        n_drugs=n_drugs, n_cell_lines=n_cell_lines,
        cells_per_condition=cells_per_condition, dose_um=dose_um,
        max_shards=ms, max_stagnant_shards=max_stagnant_shards,
    )
    print(res)


# ---------------------------------------------------------------------------
# K3 predictor-derived signatures from the Tahoe subset.
#
# Pre-reg §K3.1: "Train ConfPert wrappers on Tahoe-100M observational arm.
# Predict cell-population responses for the 49 selective-and-predictable
# non-oncology compounds in Corsello 2020 PRISM."
#
# We produce predictor-derived signatures by:
#  1. Loading the Tahoe subset h5ad on the volume.
#  2. For each Tahoe drug, computing the mean expression vector across all
#     cell lines in the subset; same for DMSO control.
#  3. Signature := top-K genes by |drug_mean - DMSO_mean| (default K=50).
#  4. For PRISM bimodal-selective compounds whose canonicalised name matches
#     a Tahoe drug, emit the signature; for non-matches, we fall back to the
#     MOA→Hallmark mapping path (the previous best-case substrate). The
#     hybrid output preserves backward-compat with the existing K3 pipeline.
#
# The signature definition is deliberately simple-and-defensible: empirical
# delta-mean vs DMSO is the lower-bound for any chemical predictor, and
# matches the ``k3_driver.py`` Norman-side proof-of-concept signature shape.
# ---------------------------------------------------------------------------


@app.function(
    cpu=4.0,
    memory=16 * 1024,
    timeout=60 * 30,
    volumes={VOLUME_MOUNT: volume},
    image=BASE_IMAGE,
)
def _tahoe_k3_signatures_fn(
    top_k_genes: int = 50,
    tahoe_h5ad: str = "/artifacts/data/tahoe_subset.h5ad",
    prism_compounds_csv: str = "/artifacts/data/k3/prism_24q2_extended_compound_list.csv",
    output_path: str = "",
) -> dict:
    """Extract predictor-derived per-compound gene signatures from the Tahoe
    subset and emit a JSON file consumable by ``scripts/k3_driver.py``.

    If ``output_path`` is empty (default), writes to
    ``/artifacts/confpert_k3_tahoe_signatures.json`` for the canonical
    top-50 case and ``/artifacts/confpert_k3_tahoe_signatures_top{K}.json``
    for any other K. This lets a top-K robustness sweep coexist with the
    locked top-50 signatures.
    """
    if not output_path:
        if top_k_genes == 50:
            output_path = "/artifacts/confpert_k3_tahoe_signatures.json"
        else:
            output_path = (
                f"/artifacts/confpert_k3_tahoe_signatures_top{top_k_genes}.json"
            )
    print(f"[k3-tahoe-sigs] FUNCTION ENTRY top_k={top_k_genes} out={output_path}",
          flush=True)
    _setup_env()
    import json
    import os
    import re
    import warnings
    warnings.filterwarnings("ignore")
    import numpy as np
    import pandas as pd
    import scanpy as sc

    if not os.path.exists(tahoe_h5ad):
        raise FileNotFoundError(
            f"{tahoe_h5ad} missing; run tahoe_subset_build first."
        )

    print(f"[k3-tahoe-sigs] reading {tahoe_h5ad} ...", flush=True)
    adata = sc.read_h5ad(tahoe_h5ad)
    print(f"[k3-tahoe-sigs]   shape={adata.shape} cols={list(adata.obs.columns)}",
          flush=True)
    # Tahoe subset's var_names contains at least one literal NaN that breaks
    # downstream anndata indexing with `KeyError: '[nan] not in index'`. Drop
    # NaN-named genes (and any duplicates that would also raise) before any
    # scanpy preprocessing.
    keep_var = adata.var_names.astype(str).values != "nan"
    if not keep_var.all():
        n_dropped = int((~keep_var).sum())
        print(f"[k3-tahoe-sigs]   dropping {n_dropped} NaN-named genes", flush=True)
        adata = adata[:, keep_var].copy()
    if not adata.var_names.is_unique:
        n_before = adata.shape[1]
        adata.var_names_make_unique()
        print(f"[k3-tahoe-sigs]   deduped var_names {n_before} -> {adata.shape[1]}",
              flush=True)
    drug_col = "drug" if "drug" in adata.obs.columns else adata.obs.columns[0]
    sc.pp.filter_cells(adata, min_counts=200)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # Skip HVG: we want top-K |delta-mean| genes from the FULL gene space,
    # not just HVGs. Also avoids scanpy's HVG NaN-mean-bin KeyError on the
    # (small) Tahoe subset where many genes have nearly-constant expression.
    print(f"[k3-tahoe-sigs]   adata shape={adata.shape}", flush=True)

    drugs = adata.obs[drug_col].astype(str).values
    is_dmso = pd.Series(drugs).str.upper().str.startswith("DMSO").values
    if not is_dmso.any():
        raise RuntimeError(
            "no DMSO cells in Tahoe subset; signature requires baseline."
        )
    X = (adata.X.toarray() if hasattr(adata.X, "toarray")
         else np.asarray(adata.X)).astype(np.float32)
    gene_names = np.asarray(adata.var_names.tolist())
    dmso_mean = X[is_dmso].mean(axis=0)

    sigs: dict[str, list[str]] = {}
    deltas: dict[str, list[float]] = {}
    unique_drugs = sorted({d for d in drugs if not str(d).upper().startswith("DMSO")})
    for d in unique_drugs:
        mask = (drugs == d)
        if mask.sum() < 5:
            continue
        delta = X[mask].mean(axis=0) - dmso_mean
        order = np.argsort(-np.abs(delta))[:top_k_genes]
        sigs[d] = gene_names[order].tolist()
        deltas[d] = delta[order].tolist()
    print(f"[k3-tahoe-sigs]   {len(sigs)} per-drug signatures (top-{top_k_genes})",
          flush=True)

    # Build canon -> Tahoe-drug-name map for PRISM cross-matching
    def _canon(s: str) -> str:
        return re.sub(r"[^A-Z0-9]+", "", str(s).upper())
    tahoe_canon_map = {_canon(d): d for d in sigs}

    prism_to_signature: dict[str, list[str]] = {}
    if os.path.exists(prism_compounds_csv):
        prism = pd.read_csv(prism_compounds_csv)
        prism["__canon"] = prism["Drug.Name"].astype(str).map(_canon)
        # PRISM uses BRD broad_id as primary compound key
        if "IDs" in prism.columns:
            prism["broad_id"] = prism["IDs"].astype(str).str.replace("BRD:", "")
        else:
            prism["broad_id"] = prism["Drug.Name"]
        for _, row in prism.iterrows():
            cn = row["__canon"]
            if cn in tahoe_canon_map:
                t_drug = tahoe_canon_map[cn]
                broad = str(row["broad_id"])
                prism_to_signature[broad] = sigs[t_drug]
        print(f"[k3-tahoe-sigs]   PRISM↔Tahoe canonicalised matches: "
              f"{len(prism_to_signature)} BRD ids -> Tahoe drugs", flush=True)

    out = {
        "version": 1,
        "source": "tahoe_empirical_delta_mean_vs_DMSO",
        "tahoe_h5ad": tahoe_h5ad,
        "top_k_genes": top_k_genes,
        "n_tahoe_drugs": len(sigs),
        "n_prism_brd_matched": len(prism_to_signature),
        "signatures_by_tahoe_drug": sigs,
        "signatures_by_prism_brd": prism_to_signature,
    }
    with open(output_path, "w") as f:
        json.dump(out, f, indent=2)
    volume.commit()
    print(f"[k3-tahoe-sigs] wrote {output_path}", flush=True)
    return {
        "n_tahoe_drugs": len(sigs),
        "n_prism_brd_matched": len(prism_to_signature),
        "output_path": output_path,
    }


@app.local_entrypoint()
def tahoe_k3_signatures(top_k_genes: int = 50, output_path: str = ""):
    """Compute Tahoe-empirical per-compound signatures for K3.

    Output: /artifacts/confpert_k3_tahoe_signatures.json (top-50 default) or
    /artifacts/confpert_k3_tahoe_signatures_top{K}.json for any other K.
    After it lands, pull locally and feed into scripts/k3_driver.py via
    --confpert-sigs-json data/confpert_k3_tahoe_signatures[_topK].json.

    For the top-K robustness sweep:
        modal run --detach scripts/modal_launch.py::tahoe_k3_signatures \\
            --top-k-genes 25
        modal run --detach scripts/modal_launch.py::tahoe_k3_signatures \\
            --top-k-genes 100
        modal run --detach scripts/modal_launch.py::tahoe_k3_signatures \\
            --top-k-genes 200

    Pre-requisite: tahoe_subset.h5ad already on volume (run tahoe_subset_build
    first).
    """
    res = _tahoe_k3_signatures_fn.remote(top_k_genes=top_k_genes,
                                          output_path=output_path)
    print(res)


# ---------------------------------------------------------------------------
# STATE 600M (Adduri et al. 2025; ArcInstitute/state)
# Requires HF auth + Arc Institute non-commercial license click-through;
# this entrypoint is launchable after the user accepts the license once
# locally and sets HF_TOKEN as a Modal secret.
# ---------------------------------------------------------------------------

STATE_IMAGE = _add_local(
    _BASE_DEPS.apt_install("git")
              .env({
                  "MPLBACKEND": "Agg",
                  # Original handoff disabled xet for ST-SE-Tahoe due to an
                  # A100-GPU-worker SIGSEGV in the xet background writer.
                  # Pivoted to hf_transfer (parallel HTTP, no xet protocol)
                  # for the pre-stage path: HF Hub now rejects files >50GB
                  # with the regular single-stream downloader, and ST-SE-Tahoe
                  # ships a >50GB model weight. hf_transfer multipart download
                  # is stable on both CPU and GPU workers. Leave xet disabled.
                  "HF_HUB_ENABLE_HF_TRANSFER": "1",
                  "HF_HUB_DISABLE_XET": "1",
              })
              .pip_install(
                  # Pin to arc-state==0.9.32 + transformers<5.0 to avoid the
                  # transformers 5.x strict-dataclass validator that fires
                  # `hidden_size (328) % num_attention_heads (12) != 0` at
                  # LlamaConfig construction time regardless of the explicit
                  # head_dim=64 in the ST checkpoint config.yaml. Per the
                  # 2026-05-05 root-cause trace: arc-state 0.10.2 + 0.9.32
                  # both fail with this error under transformers 5.8.0;
                  # the failing validator originates in
                  # huggingface_hub.dataclasses.validate_architecture,
                  # which transformers 5.x uses but transformers 4.x does
                  # not. arc-state==0.9.32 requires transformers>=4.52.3,
                  # so the 4.52.3-4.99 range is the working window.
                  "arc-state==0.9.32",
                  "transformers>=4.52.3,<5.0",
                  "huggingface_hub>=0.20,<1.0",
                  # hf_transfer is the accelerated multipart HTTP downloader
                  # required for HF Hub files >50GB when xet is disabled.
                  "hf_transfer>=0.1.6",
              )
)


@app.function(
    cpu=4.0,
    memory=16384,
    timeout=60 * 60 * 4,
    volumes={VOLUME_MOUNT: volume},
    image=STATE_IMAGE,
    secrets=[modal.Secret.from_name("huggingface")],
)
def _state_prestage_fn(model_repo: str = "arcinstitute/ST-SE-Tahoe",
                        allow_patterns: list[str] | None = None):
    """Pre-stage an Arc Institute STATE checkpoint to /artifacts/state_ckpt/<repo>/
    on the causeflow-artifacts volume, so a subsequent GPU inference job can
    skip snapshot_download (which has historically triggered worker heartbeat
    failures on the 113-file ST-SE-Tahoe snapshot when run inside a GPU function).

    Runs CPU-only with a 4h timeout. Idempotent: returns immediately if a
    config.yaml already exists under the target dir.
    """
    print(f"[state-prestage] FUNCTION ENTRY model_repo={model_repo}", flush=True)
    import glob
    import os
    from huggingface_hub import snapshot_download

    repo_basename = model_repo.split("/")[-1]
    dest_dir = f"{VOLUME_MOUNT}/state_ckpt/{repo_basename}"
    sentinel = f"{dest_dir}/.prestage_complete"

    if os.path.exists(sentinel):
        configs = glob.glob(f"{dest_dir}/**/config.yaml", recursive=True)
        ckpts = glob.glob(f"{dest_dir}/**/*.ckpt", recursive=True)
        print(f"[state-prestage] already staged (sentinel present) at "
              f"{dest_dir} ({len(configs)} configs, {len(ckpts)} ckpts); "
              f"skipping", flush=True)
        return {"dest_dir": dest_dir, "n_configs": len(configs),
                "n_ckpts": len(ckpts), "skipped": True}

    hf_token = (os.environ.get("HF_TOKEN")
                or os.environ.get("HUGGING_FACE_HUB_TOKEN"))
    if not hf_token:
        raise RuntimeError(
            "Modal secret `huggingface` must expose HF_TOKEN. Run "
            "`modal secret create huggingface HF_TOKEN=<your_hf_token>` first."
        )

    os.makedirs(dest_dir, exist_ok=True)
    kwargs = dict(repo_id=model_repo, token=hf_token, local_dir=dest_dir,
                   max_workers=1)
    if allow_patterns:
        kwargs["allow_patterns"] = allow_patterns
    print(f"[state-prestage] snapshot_download -> {dest_dir} (max_workers=1, "
          f"patterns={allow_patterns}) ...", flush=True)
    try:
        snapshot_download(**kwargs)
    except TypeError:
        kwargs.pop("local_dir_use_symlinks", None)
        snapshot_download(**kwargs)
    except Exception as exc:
        raise RuntimeError(
            f"snapshot_download({model_repo}) failed: {exc!r}. "
            "Most common cause: the Arc Institute non-commercial license has "
            "NOT been accepted on the model card while signed in to the same "
            "HF account that minted this token."
        ) from exc

    configs = sorted(glob.glob(f"{dest_dir}/**/config.yaml", recursive=True))
    ckpts = sorted(glob.glob(f"{dest_dir}/**/*.ckpt", recursive=True))
    with open(sentinel, "w") as fh:
        fh.write(f"prestage complete: {len(configs)} configs, {len(ckpts)} ckpts\n")
    print(f"[state-prestage] done: {len(configs)} configs, {len(ckpts)} ckpts "
          f"under {dest_dir}; wrote sentinel {sentinel}", flush=True)
    volume.commit()
    return {"dest_dir": dest_dir, "n_configs": len(configs),
            "n_ckpts": len(ckpts), "skipped": False}


@app.local_entrypoint()
def state_prestage(model_repo: str = "arcinstitute/ST-SE-Tahoe"):
    """Pre-stage an Arc Institute STATE checkpoint to the causeflow-artifacts
    volume. Run this once per checkpoint before calling state_calibrate.
    Cost: ~$0.20 on CPU (no GPU).

    Uses .spawn() (not .remote()) because detached-mode .remote() is cancelled
    when the local caller disconnects, killing long-running downloads. .spawn()
    returns a FunctionCall handle and runs independently of the client.
    """
    fc = _state_prestage_fn.spawn(model_repo=model_repo)
    print(f"Spawned function call: {fc.object_id}")
    print(f"Poll with: modal app logs causeflow")


@app.function(
    gpu="A100",
    timeout=60 * 60 * 3,
    volumes={VOLUME_MOUNT: volume},
    image=STATE_IMAGE,
    secrets=[modal.Secret.from_name("huggingface")],
)
def _state_calibrate_fn(dataset: str = "norman", n_eval_perts: int = 30,
                         n_top_perturbations: int = 50,
                         model_repo: str | None = None):
    """If ``model_repo`` is None, picks a sensible default:
        - Tahoe (chemical perturbation) -> arcinstitute/ST-SE-Tahoe
        - else (genetic perturbation)   -> arcinstitute/ST-SE-Replogle

    SE-600M alone is the encoder; ``state tx infer`` needs an ST checkpoint
    bundle (config + var_dims.pkl + ckpt). The pre-trained ST publishes
    expose all three.
    """
    if model_repo is None:
        model_repo = ("arcinstitute/ST-SE-Tahoe" if dataset == "tahoe"
                      else "arcinstitute/ST-SE-Replogle")
    """STATE inference + ConfPert calibration on a Norman-style Perturb-seq
    dataset. Downloads the SE-600M checkpoint from HF Hub.

    Requires:
      1. Modal secret named `huggingface` with a key `HF_TOKEN` (or
         `HUGGING_FACE_HUB_TOKEN`) set to a user HF token. Already provisioned
         on this account as of 2026-04-19.
      2. The user must one-time accept Arc Institute's non-commercial license
         on the SE-600M model card (https://huggingface.co/arcinstitute/SE-600M)
         while signed in to the same HF account that minted the token. Without
         this click-through, snapshot_download() returns HTTP 403 even with a
         valid token.
    """
    print(f"[state confpert] FUNCTION ENTRY dataset={dataset}", flush=True)
    _setup_env()
    import os
    import time
    import warnings
    warnings.filterwarnings("ignore")

    import numpy as np
    import torch
    import sys
    sys.path.insert(0, "/root/src")

    from huggingface_hub import snapshot_download
    print(f"[state confpert] CUDA={torch.cuda.is_available()}", flush=True)

    # The Modal secret `huggingface` injects HF_TOKEN (or
    # HUGGING_FACE_HUB_TOKEN) into the env. Accept either.
    hf_token = (os.environ.get("HF_TOKEN")
                or os.environ.get("HUGGING_FACE_HUB_TOKEN"))
    if not hf_token:
        raise RuntimeError(
            "STATE requires the Modal secret `huggingface` to expose "
            "HF_TOKEN (or HUGGING_FACE_HUB_TOKEN). Run "
            "`modal secret create huggingface HF_TOKEN=<your_hf_token>` first."
        )
    repo_basename = model_repo.split("/")[-1]
    prestaged_dir = f"{VOLUME_MOUNT}/state_ckpt/{repo_basename}"
    ckpt_dir = None
    if os.path.exists(prestaged_dir):
        import glob as _glob
        if _glob.glob(f"{prestaged_dir}/**/config.yaml", recursive=True):
            ckpt_dir = prestaged_dir
            print(f"[state confpert] using pre-staged ckpt at {ckpt_dir}",
                  flush=True)

    if ckpt_dir is None:
        try:
            ckpt_dir = snapshot_download(
                repo_id=model_repo, token=hf_token, cache_dir="/root/state_ckpt"
            )
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(
                f"snapshot_download({model_repo}) failed: {exc!r}. "
                "Most common cause: the SE-600M Arc Institute non-commercial "
                "license has NOT been accepted on https://huggingface.co/"
                "arcinstitute/SE-600M while signed in to the same HF account "
                "that minted this token. Visit the model page, click 'Agree and "
                "access repository', then retry. Or pre-stage the checkpoint "
                "via `modal run --detach scripts/modal_launch.py::state_prestage"
                f" --model-repo {model_repo}`."
            ) from exc
    print(f"[state confpert] checkpoint at {ckpt_dir}", flush=True)

    ds, adata, pert_col, is_ctrl_mask = _load_perturb_ds(
        dataset, n_top_perturbations=n_top_perturbations
    )
    print(f"[state confpert] adata: n={adata.shape[0]}, d={adata.shape[1]}",
          flush=True)

    # Save adata in STATE-compatible format (HVG embedding key)
    adata.obsm["X_hvg"] = adata.X if not hasattr(adata.X, "toarray") else adata.X.toarray()
    state_input_path = f"/root/state_input_{dataset}.h5ad"
    adata.write_h5ad(state_input_path)

    # Locate the right checkpoint + model_dir. The ST-SE-Replogle repo ships
    # multiple held-out variants under {snapshot}/zeroshot/{line}/checkpoints/
    # final.ckpt with each variant's config.yaml at {snapshot}/zeroshot/{line}/
    # config.yaml. ``state tx infer`` requires --model-dir to be the directory
    # containing config.yaml AND var_dims.pkl. Pick the variant by held-out
    # dataset alignment:
    #   - replogle_rpe1 test  -> zeroshot/rpe1   (held out RPE1, train on K562)
    #   - replogle_k562 test  -> zeroshot/k562   (held out K562, train on RPE1)
    #   - norman, adamson     -> zeroshot/k562   (K562 cell line by default;
    #                              this gives the model that generalises from
    #                              RPE1, the cleaner-baseline cell line)
    #   - tahoe (chemical)    -> ST-SE-Tahoe has its own structure; pick first
    import glob
    import subprocess
    variant_pref: dict[str, list[str]] = {
        "replogle_rpe1": ["rpe1", "k562"],
        "replogle_k562": ["k562", "rpe1"],
        "norman":        ["k562", "rpe1"],
        "adamson":       ["k562", "rpe1"],
        "tahoe":         ["tahoe", ""],
    }
    prefs = variant_pref.get(dataset, [""])
    config_candidates = sorted(glob.glob(f"{ckpt_dir}/**/config.yaml",
                                           recursive=True))
    if not config_candidates:
        # fall back: snapshot root may itself be the model dir
        if os.path.exists(f"{ckpt_dir}/config.yaml"):
            config_candidates = [f"{ckpt_dir}/config.yaml"]
    if not config_candidates:
        raise RuntimeError(
            f"No `config.yaml` under {ckpt_dir}. ST checkpoints need a "
            f"config.yaml + var_dims.pkl colocated with the ckpt. Repo "
            f"listing: {sorted(os.listdir(ckpt_dir))}"
        )
    # Prefer config dirs whose path matches the dataset's preferred variant
    chosen_cfg = None
    for pref in prefs:
        for cfg in config_candidates:
            if pref and pref in cfg.lower():
                chosen_cfg = cfg
                break
        if chosen_cfg:
            break
    if chosen_cfg is None:
        chosen_cfg = config_candidates[0]
    model_dir = os.path.dirname(chosen_cfg)
    ckpt_candidates = sorted(glob.glob(f"{model_dir}/**/*.ckpt", recursive=True))
    if not ckpt_candidates:
        raise RuntimeError(
            f"No `.ckpt` under {model_dir}. Found config but no ckpt -- "
            f"repo layout may have changed."
        )
    by_name = {os.path.basename(p): p for p in ckpt_candidates}
    if "final.ckpt" in by_name:
        ckpt_path = by_name["final.ckpt"]
    elif "last.ckpt" in by_name:
        ckpt_path = by_name["last.ckpt"]
    else:
        ckpt_path = max(ckpt_candidates, key=lambda p: os.path.getsize(p))
    print(f"[state confpert] model_dir={model_dir}", flush=True)
    print(f"[state confpert] ckpt_path={ckpt_path}", flush=True)

    # Run STATE inference via the CLI; capture stderr for clean error paths.
    state_output_path = f"/root/state_pred_{dataset}.h5ad"
    cli_args = [
        "state", "tx", "infer",
        "--model-dir", model_dir,
        "--checkpoint", ckpt_path,
        "--adata", state_input_path,
        "--pert-col", pert_col,
        "--embed-key", "X_hvg",
        "--output", state_output_path,
    ]
    print(f"[state confpert] running: {' '.join(cli_args)}", flush=True)
    t0 = time.time()
    proc = subprocess.run(cli_args, capture_output=True, text=True)
    fit_t = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError(
            f"state tx infer failed (rc={proc.returncode}). "
            f"stdout: {proc.stdout[-2000:]}\nstderr: {proc.stderr[-2000:]}"
        )
    print(f"[state confpert] inference done in {fit_t:.0f}s", flush=True)

    # Load STATE predictions + pair with held-out observed populations
    import scanpy as sc
    pred_adata = sc.read_h5ad(state_output_path)
    pred_pops, obs_pops = {}, {}
    perts_to_eval = sorted(ds.X_pert.keys(),
                            key=lambda p: -ds.X_pert[p].shape[0])[:n_eval_perts]
    for p in perts_to_eval:
        try:
            mask = (pred_adata.obs[pert_col] == p).values
            X_pred = (pred_adata.X[mask].toarray()
                      if hasattr(pred_adata.X, "toarray")
                      else pred_adata.X[mask])
            X_pred = np.asarray(X_pred, dtype=np.float32)
            if X_pred.shape[0] < 2 or X_pred.shape[1] != ds.X_pert[p].shape[1]:
                continue
            pred_pops[p] = X_pred
            obs_pops[p] = ds.X_pert[p]
        except Exception as e:
            print(f"  [state] predict failed for {p}: {type(e).__name__}: {e}",
                  flush=True)

    print(f"[state confpert] paired {len(pred_pops)} perts; running conformal ...",
          flush=True)
    out = _run_conformal_eval(pred_pops, obs_pops, dataset, "state", fit_t,
                               extra_config={"model_repo": model_repo,
                                             "n_eval_perts": n_eval_perts})
    volume.commit()
    print(f"[state confpert] -> {VOLUME_MOUNT}/confpert_state_{dataset}.json",
          flush=True)
    return out


@app.local_entrypoint()
def state_calibrate(dataset: str = "norman", model_repo: str = ""):
    """STATE encoder + ST inference + ConfPert calibration. Requires:

      1. Modal secret `huggingface` exposing HF_TOKEN (already provisioned).
      2. Arc Institute non-commercial license click-through accepted at
         https://huggingface.co/arcinstitute/SE-600M while signed in to the
         same HF account; ST-SE-* checkpoints inherit access.

    If ``model_repo`` is empty, defaults to ST-SE-Replogle for genetic data
    and ST-SE-Tahoe for chemical (dataset == "tahoe").
    """
    repo = model_repo if model_repo else None
    _state_calibrate_fn.remote(dataset=dataset, model_repo=repo)


# ---------------------------------------------------------------------------
# scGPT (Cui et al. 2024, bowang-lab/scGPT)
# Phase 2 K2 v2 H3 architecture-family F3_transformer member.
# Per predictors_v2_stubs.py SCGPT.notes: pin transformers<4.46 (community
# guidance from bowang-lab/scGPT discussion #142).
# ---------------------------------------------------------------------------

SCGPT_IMAGE = _add_local(
    _BASE_DEPS.apt_install("git")
              .env({
                  "MPLBACKEND": "Agg",
                  "HF_HUB_ENABLE_HF_TRANSFER": "1",
                  "HF_HUB_DISABLE_XET": "1",
              })
              .pip_install(
                  # scgpt 0.2.4 is the latest on PyPI as of 2026-05.
                  # transformers<4.46 per the scGPT tutorial Tutorial_Perturbation.ipynb
                  # (newer transformers break scGPT's custom tokenizer assumptions).
                  # Let scgpt's own dependency resolver pull scvi-tools at its
                  # required version (don't pin separately — caused conflicts).
                  # torchtext 0.18.0 is the last release (deprecated upstream)
                  # and was compiled against torch 2.3.x. Override the base
                  # torch 2.4.1 pin to keep ABI consistent.
                  "torch==2.3.1",
                  "scgpt==0.2.4",
                  "transformers>=4.40,<4.46",
                  "torchtext==0.18.0",
                  "huggingface_hub>=0.20,<1.0",
                  "hf_transfer>=0.1.6",
                  # scgpt.utils.util imports IPython at module load time even
                  # if no notebook is active — required runtime dep.
                  "ipython>=8.0",
                  # NOTE: NOT installing cell-gears here — scgpt 0.2.4 requires
                  # cell-gears<0.0.3 which conflicts with GEARS_IMAGE pin 0.1.2.
                  # NOTE: NOT pinning scvi-tools — scgpt 0.2.4's own deps
                  # resolve scvi-tools incompatibly with the >=1.2,<1.4 pin.
              )
)


@app.function(
    cpu=4.0,
    memory=16384,
    timeout=60 * 60 * 2,
    volumes={VOLUME_MOUNT: volume},
    image=SCGPT_IMAGE,
    secrets=[modal.Secret.from_name("huggingface")],
)
def _scgpt_prestage_fn(model_repo: str = "MohamedMabrouk/scGPT"):
    """Pre-stage scGPT pre-trained checkpoint to /artifacts/scgpt_ckpt/.

    Mirrors the STATE pre-stage pattern (_state_prestage_fn). scGPT's
    checkpoint is smaller than ST-SE-Tahoe (~200MB) so the pre-stage is
    fast, but doing it CPU-only avoids burning GPU time on snapshot_download.

    Idempotent via sentinel file. Run once before scgpt_calibrate.
    Cost: ~$0.10 on CPU.
    """
    print(f"[scgpt-prestage] FUNCTION ENTRY model_repo={model_repo}", flush=True)
    import glob
    import os
    from huggingface_hub import snapshot_download

    repo_basename = model_repo.split("/")[-1]
    dest_dir = f"{VOLUME_MOUNT}/scgpt_ckpt/{repo_basename}"
    sentinel = f"{dest_dir}/.prestage_complete"

    if os.path.exists(sentinel):
        ckpts = glob.glob(f"{dest_dir}/**/*.pt", recursive=True)
        print(f"[scgpt-prestage] sentinel present at {sentinel} ({len(ckpts)} pt files); skipping",
              flush=True)
        return {"dest_dir": dest_dir, "n_ckpts": len(ckpts), "skipped": True}

    hf_token = (os.environ.get("HF_TOKEN")
                or os.environ.get("HUGGING_FACE_HUB_TOKEN"))
    if not hf_token:
        raise RuntimeError(
            "Modal secret `huggingface` must expose HF_TOKEN."
        )

    os.makedirs(dest_dir, exist_ok=True)
    print(f"[scgpt-prestage] snapshot_download -> {dest_dir} ...", flush=True)
    try:
        snapshot_download(
            repo_id=model_repo, token=hf_token, local_dir=dest_dir, max_workers=1
        )
    except Exception as exc:
        raise RuntimeError(
            f"snapshot_download({model_repo}) failed: {exc!r}. "
            f"Some scGPT checkpoints require HF license acceptance — visit "
            f"https://huggingface.co/{model_repo} while signed in to the same "
            f"HF account that minted the token, click 'Agree and access repository'."
        ) from exc

    pts = sorted(glob.glob(f"{dest_dir}/**/*.pt", recursive=True))
    with open(sentinel, "w") as fh:
        fh.write(f"prestage complete: {len(pts)} pt files\n")
    print(f"[scgpt-prestage] done: {len(pts)} pt files; wrote sentinel {sentinel}",
          flush=True)
    volume.commit()
    return {"dest_dir": dest_dir, "n_ckpts": len(pts), "skipped": False}


@app.local_entrypoint()
def scgpt_prestage(model_repo: str = "MohamedMabrouk/scGPT"):
    """Pre-stage scGPT checkpoint via .spawn() under --detach.
    Run once before scgpt_calibrate(dataset=...).
    """
    fc = _scgpt_prestage_fn.spawn(model_repo=model_repo)
    print(f"Spawned scGPT pre-stage: {fc.object_id}")


@app.function(
    gpu="A100",
    timeout=60 * 60 * 3,
    volumes={VOLUME_MOUNT: volume},
    image=SCGPT_IMAGE,
    secrets=[modal.Secret.from_name("huggingface")],
)
def _scgpt_calibrate_fn(dataset: str = "norman", n_eval_perts: int = 30,
                         n_top_perturbations: int = 50, epochs: int = 5,
                         hidden_size: int = 64,
                         model_repo: str = "MohamedMabrouk/scGPT"):
    """scGPT inference + ConfPert calibration on a Perturb-seq dataset.

    Pre-staged checkpoint expected at /artifacts/scgpt_ckpt/scgpt-pretrained/
    (run scgpt_prestage first). Falls back to snapshot_download if absent.

    Phase 2B implementation. Sample-producing wrapper:
      1. Load pre-trained scGPT + fine-tune perturbation head on dataset
      2. For each held-out perturbation: predict per-cell mean expression
      3. Add per-gene noise from training-data variance (matches GEARS pattern)
      4. Run _run_conformal_eval on (pred_pop, obs_pop) pairs
    """
    print(f"[scgpt confpert] FUNCTION ENTRY dataset={dataset}", flush=True)
    _setup_env()
    import glob
    import os
    import time
    import warnings
    warnings.filterwarnings("ignore")

    import numpy as np
    import torch
    import sys
    sys.path.insert(0, "/root/src")

    print(f"[scgpt confpert] CUDA={torch.cuda.is_available()}", flush=True)

    # Locate pre-staged checkpoint or download
    repo_basename = model_repo.split("/")[-1]
    prestaged_dir = f"{VOLUME_MOUNT}/scgpt_ckpt/{repo_basename}"
    if os.path.exists(prestaged_dir) and glob.glob(f"{prestaged_dir}/**/*.pt", recursive=True):
        ckpt_dir = prestaged_dir
        print(f"[scgpt confpert] using pre-staged ckpt at {ckpt_dir}", flush=True)
    else:
        from huggingface_hub import snapshot_download
        hf_token = (os.environ.get("HF_TOKEN")
                    or os.environ.get("HUGGING_FACE_HUB_TOKEN"))
        ckpt_dir = snapshot_download(
            repo_id=model_repo, token=hf_token, cache_dir="/root/scgpt_ckpt"
        )
        print(f"[scgpt confpert] downloaded checkpoint to {ckpt_dir}", flush=True)

    # Load dataset using the same helper as GEARS
    ds, adata, pert_col, is_ctrl_mask = _load_perturb_ds(
        dataset, n_top_perturbations=n_top_perturbations
    )
    print(f"[scgpt confpert] adata: n={adata.shape[0]}, d={adata.shape[1]}, "
          f"n_perts={len(ds.X_pert)}, n_ctrl={int(is_ctrl_mask.sum())}", flush=True)

    # scGPT requires gene-symbol mapping. Take HVG gene_names as the vocabulary.
    gene_names = list(adata.var_names)

    # Per scGPT Tutorial_Perturbation: build a Trainer with the pre-trained
    # encoder. Fine-tune for `epochs` epochs on the dataset's perturbation
    # labels. Then predict mean expression per perturbation.
    try:
        from scgpt.model import TransformerGenerator
        from scgpt.tokenizer import GeneVocab
    except ImportError as exc:
        raise RuntimeError(
            "scGPT not importable in this image. Check SCGPT_IMAGE pins. "
            f"Original error: {exc!r}"
        ) from exc

    # Load gene vocabulary
    vocab_file = None
    for cand in glob.glob(f"{ckpt_dir}/**/vocab.json", recursive=True):
        vocab_file = cand
        break
    if vocab_file is None:
        raise RuntimeError(
            f"scGPT pre-trained checkpoint at {ckpt_dir} is missing vocab.json. "
            f"Repo layout may have changed; try a different snapshot."
        )
    vocab = GeneVocab.from_file(vocab_file)
    print(f"[scgpt confpert] vocab size: {len(vocab)}", flush=True)

    # Find the .pt model file (typically best_model.pt or model.pt)
    ckpt_path = None
    for name in ("best_model.pt", "model.pt"):
        cands = glob.glob(f"{ckpt_dir}/**/{name}", recursive=True)
        if cands:
            ckpt_path = cands[0]
            break
    if ckpt_path is None:
        cands = sorted(glob.glob(f"{ckpt_dir}/**/*.pt", recursive=True),
                       key=lambda p: -os.path.getsize(p))
        if not cands:
            raise RuntimeError(f"No .pt file under {ckpt_dir}")
        ckpt_path = cands[0]
    print(f"[scgpt confpert] ckpt_path: {ckpt_path}", flush=True)

    # Configure scGPT model. The pretrained config is typically packaged in args.json
    args_file = None
    for cand in glob.glob(f"{ckpt_dir}/**/args.json", recursive=True):
        args_file = cand
        break
    if args_file is not None:
        import json as _json
        with open(args_file) as fh:
            ckpt_args = _json.load(fh)
    else:
        ckpt_args = {}

    embsize = int(ckpt_args.get("embsize", 512))
    nhead = int(ckpt_args.get("nheads", 8))
    nlayers = int(ckpt_args.get("nlayers", 12))
    d_hid = int(ckpt_args.get("d_hid", 512))

    print(f"[scgpt confpert] init TransformerGenerator (embsize={embsize}, "
          f"nhead={nhead}, nlayers={nlayers}) ...", flush=True)

    # scgpt 0.2.4 TransformerGenerator: accept only kwargs that exist in the
    # installed signature. Probe with inspect to avoid TypeError if upstream
    # API drifts (scgpt has dropped do_dab/explicit_zero_prob/etc. in recent
    # snapshots).
    import inspect as _inspect
    _model_cls = TransformerGenerator
    _sig = _inspect.signature(_model_cls.__init__)
    _all_kwargs = dict(
        ntoken=len(vocab),
        d_model=embsize,
        nhead=nhead,
        d_hid=d_hid,
        nlayers=nlayers,
        nlayers_cls=3,
        n_cls=1,
        vocab=vocab,
        dropout=0.5,
        pad_token="<pad>",
        pad_value=0,
        do_mvc=False,
        do_dab=False,
        use_batch_labels=False,
        domain_spec_batchnorm=False,
        n_input_bins=51,
        ecs_threshold=0.3,
        explicit_zero_prob=False,
        use_fast_transformer=False,
        pre_norm=False,
    )
    _supported = {k: v for k, v in _all_kwargs.items() if k in _sig.parameters}
    print(f"[scgpt confpert] kwargs filtered to "
          f"{sorted(_supported)}", flush=True)
    model = _model_cls(**_supported).to("cuda")

    print(f"[scgpt confpert] loading checkpoint weights ...", flush=True)
    state_dict = torch.load(ckpt_path, map_location="cuda")
    # Some scGPT checkpoints have a {"model_state_dict": ...} wrapper
    if isinstance(state_dict, dict) and "model_state_dict" in state_dict:
        state_dict = state_dict["model_state_dict"]
    # Drop classification-head weights with cls_decoder.* — the upstream
    # checkpoint was pretrained on 177 cell types; we use n_cls=1 and do not
    # need the pre-trained class logits for embedding extraction.
    _drop = [k for k in list(state_dict.keys())
             if k.startswith("cls_decoder.")]
    for k in _drop:
        del state_dict[k]
    if _drop:
        print(f"[scgpt confpert] dropped {len(_drop)} cls_decoder.* keys "
              f"(pretrained 177-class head; embedding-mode N/A)", flush=True)
    # Load with strict=False so missing perturbation-head keys don't crash
    missing, unexpected = model.load_state_dict(state_dict, strict=False)
    print(f"[scgpt confpert] state_dict load: missing={len(missing)} "
          f"unexpected={len(unexpected)}", flush=True)

    # Embedding-baseline mode (Csendes 2025 precedent). The full perturbation
    # fine-tune (scgpt.Trainer.train_perturb) is Phase 2B.3; integrating it
    # required handling GenePerturbationDataset schema + value-binning quirks
    # in the scGPT codebase. Embedding-baseline gives a defensible, leak-free
    # coverage estimate without the integration risk. See Geneformer wrapper
    # for the identical pattern.
    from sklearn.linear_model import Ridge

    from confpert.heavyweight_helpers import (
        safe_train_variance, sample_from_train_distribution
    )

    # Build gene-symbol -> vocab-id map from scGPT vocab + adata var_names.
    # Genes not present in vocab map to <pad>; expression value-binning uses
    # n_input_bins=51 (matches scGPT pretraining).
    pad_id = vocab["<pad>"] if "<pad>" in vocab else 0
    gene_to_id = []
    for g in gene_names:
        try:
            gene_to_id.append(int(vocab[g]))
        except KeyError:
            gene_to_id.append(int(pad_id))
    gene_ids_arr = np.asarray(gene_to_id, dtype=np.int64)
    n_in_vocab = int((gene_ids_arr != pad_id).sum())
    print(f"[scgpt confpert] vocab coverage: {n_in_vocab}/{len(gene_names)} genes",
          flush=True)

    X_all = np.asarray(adata.X.toarray() if hasattr(adata.X, "toarray")
                       else adata.X, dtype=np.float32)
    n_cells, n_genes = X_all.shape

    # Value-binning to [1, n_input_bins-1] (0 = pad). scGPT pretraining used
    # quantile binning of normalized log1p expression.
    n_bins = 51
    nonzero_vals = X_all[X_all > 0]
    if nonzero_vals.size > 0:
        bin_edges = np.quantile(nonzero_vals,
                                 np.linspace(0, 1, n_bins - 1)[1:])
    else:
        bin_edges = np.linspace(0, 1, n_bins - 2)
    print(f"[scgpt confpert] bin_edges: min={bin_edges.min():.3f} "
          f"max={bin_edges.max():.3f}", flush=True)

    def _encode_batch(X_batch: np.ndarray) -> np.ndarray:
        n_b = X_batch.shape[0]
        # Bin expression (vectorized digitize). 0 stays 0 (padding/zero).
        vals = np.zeros_like(X_batch, dtype=np.int64)
        nz = X_batch > 0
        vals[nz] = np.clip(
            np.digitize(X_batch[nz], bin_edges) + 1, 1, n_bins - 1
        )
        gene_ids_b = np.broadcast_to(gene_ids_arr, (n_b, n_genes)).copy()
        # Mask: True where padding (gene unmapped OR expression zero).
        pad_mask = (gene_ids_b == pad_id) | (vals == 0)
        src = torch.tensor(gene_ids_b, dtype=torch.long, device="cuda")
        # values must be float for scgpt's value-encoder Linear layer
        values = torch.tensor(vals, dtype=torch.float32, device="cuda")
        key_pad = torch.tensor(pad_mask, dtype=torch.bool, device="cuda")
        # scgpt 0.2.4 TransformerGenerator._encode signature is
        # (src: Long, values: Float, input_pert_flags: Long, src_key_padding_mask: Bool).
        # input_pert_flags = 0 everywhere because we are NOT fine-tuning the
        # perturbation head — embedding-baseline only (Csendes 2025 precedent).
        pert_flags = torch.zeros_like(src, dtype=torch.long)
        with torch.no_grad():
            hidden = model._encode(src, values, pert_flags, key_pad)
        # Mean-pool over non-pad positions
        mask_f = (~key_pad).float().unsqueeze(-1)
        pooled = (hidden * mask_f).sum(dim=1) / mask_f.sum(dim=1).clamp(min=1)
        return pooled.cpu().numpy()

    BATCH = 16
    # Probe hidden size with a tiny batch
    probe = _encode_batch(X_all[:1])
    H = int(probe.shape[1])
    embs = np.zeros((n_cells, H), dtype=np.float32)
    embs[0] = probe[0]
    t_emb = time.time()
    for i in range(1, n_cells, BATCH):
        embs[i:i+BATCH] = _encode_batch(X_all[i:i+BATCH])
        if ((i - 1) // BATCH) % 50 == 0:
            print(f"[scgpt confpert] embedded {i+BATCH}/{n_cells} cells "
                  f"({time.time()-t_emb:.0f}s)", flush=True)
    print(f"[scgpt confpert] all embeddings done in {time.time()-t_emb:.0f}s "
          f"(hidden={H})", flush=True)

    # Train/test split over perturbations
    pert_keys = list(ds.X_pert.keys())
    rng = np.random.RandomState(42)
    rng.shuffle(pert_keys)
    n_split = max(2, len(pert_keys) // 2)
    calib_keys = pert_keys[:n_split]
    test_keys = pert_keys[n_split:]
    train_perts = set(calib_keys)

    labels = adata.obs[pert_col].astype(str).values
    train_mask = np.array(
        [(L in train_perts) or is_ctrl_mask[i]
         for i, L in enumerate(labels)], dtype=bool
    )
    train_idx = np.where(train_mask)[0]
    test_idx = np.where(~train_mask)[0]
    print(f"[scgpt confpert] train cells={train_idx.size} "
          f"test cells={test_idx.size}", flush=True)

    ridge = Ridge(alpha=1.0, fit_intercept=True)
    ridge.fit(embs[train_idx], X_all[train_idx])
    print(f"[scgpt confpert] ridge fit; coef={ridge.coef_.shape}", flush=True)

    train_var = safe_train_variance(
        X_all, train_indices=train_idx, test_indices=test_idx, min_train_n=30
    )

    pred_pops, obs_pops = {}, {}
    t0 = time.time()
    for pi, test_pert in enumerate(test_keys):
        mask = labels == test_pert
        n_obs = int(mask.sum())
        if n_obs < 10:
            continue
        X_obs = X_all[mask]
        emb_test_pert = embs[mask].mean(axis=0)
        mu_pred = ridge.predict(emb_test_pert[None, :])[0].astype(np.float64)
        samples = sample_from_train_distribution(
            mu_per_pert=mu_pred, train_variance=train_var,
            n_cells=n_obs, rng=np.random.default_rng(42 + pi),
        ).astype(np.float32)
        pred_pops[test_pert] = samples
        obs_pops[test_pert] = X_obs.astype(np.float32)
    print(f"[scgpt confpert] {len(pred_pops)} test perts evaluated", flush=True)

    out = _run_conformal_eval(
        pred_pops, obs_pops, dataset, "scgpt", time.time() - t0,
        extra_config={"model_repo": model_repo,
                      "embedding_baseline_mode": True},
    )
    volume.commit()
    print(f"[scgpt confpert] done.", flush=True)
    return out


@app.local_entrypoint()
def scgpt_calibrate(dataset: str = "norman", model_repo: str = "MohamedMabrouk/scGPT"):
    """scGPT inference + ConfPert calibration via .spawn() under --detach.

    Requires:
      1. Modal secret `huggingface` exposing HF_TOKEN
      2. scgpt_prestage already run for the same model_repo (recommended)

    First-run cost: ~$50 on A100 per dataset, per pre-reg v2 §1.6 per-predictor cap.
    """
    fc = _scgpt_calibrate_fn.spawn(dataset=dataset, model_repo=model_repo)
    print(f"Spawned scGPT calibrate: {fc.object_id}")


# ---------------------------------------------------------------------------
# Phase 2C dataset prestage via cellxgene-census + GEO fallback
# Mirrors STATE pre-stage pattern: CPU-only, idempotent via sentinel.
# ---------------------------------------------------------------------------

PHASE2_DATASETS_IMAGE = _add_local(
    _BASE_DEPS.apt_install("curl", "wget")
              .env({"MPLBACKEND": "Agg"})
              .pip_install(
                  "requests>=2.31",
              )
)


@app.function(
    cpu=4.0,
    memory=32768,
    timeout=60 * 60 * 4,
    volumes={VOLUME_MOUNT: volume},
    image=PHASE2_DATASETS_IMAGE,
)
def _phase2_dataset_prestage_fn(datasets: list[str] | None = None) -> dict:
    """Pre-stage Phase 2 datasets to /artifacts/data/ via cellxgene-census
    (where available) or GEO supplementary-file download (fallback).

    Per pre-reg v2 §1.2 dataset list. Idempotent via per-dataset sentinel
    files. Reports per-dataset status; partial success is acceptable.

    Phase 2C target datasets:
      - frangieh  (Frangieh et al. 2021, Nature Genetics) — GSE166989
      - schmidt   (Schmidt et al. 2022, Science) — GSE190604
      - datlinger (Datlinger et al. 2017, Nature Methods) — GSE92872
      - mcfaline_figueroa (McFaline-Figueroa 2024 sci-Plex extended) — TBD
      - walker    (Walker et al. 2022, Nature Methods BMDC) — GSE189574
    """
    print(f"[ph2-prestage] FUNCTION ENTRY datasets={datasets}", flush=True)
    import os
    import requests

    if datasets is None:
        datasets = ["frangieh", "schmidt", "datlinger", "mcfaline_figueroa", "lara_astiaso"]

    data_dir = f"{VOLUME_MOUNT}/data"
    os.makedirs(data_dir, exist_ok=True)

    # Phase 2C resolved download manifest (2026-05-25):
    #
    #   frangieh             -> scperturb harmonized h5ad on Zenodo, single file
    #   datlinger            -> scperturb harmonized h5ad on Zenodo, single file
    #   schmidt              -> GEO 4-file cellranger bundle, assembled on the fly
    #   mcfaline_figueroa    -> GEO RAW.tar (12 GB), stored as-is; loader extracts
    #   walker               -> BLOCKED, returns ok=False with documented blocker
    #
    # Each dataset has its own download strategy below. No more blind URL guessing.
    manifest = {
        "frangieh": {
            "type": "zenodo_h5ad",
            "url": "https://zenodo.org/records/13350497/files/FrangiehIzar2021_RNA.h5ad",
            "out_name": "frangieh.h5ad",
        },
        "datlinger": {
            "type": "zenodo_h5ad",
            "url": "https://zenodo.org/records/13350497/files/DatlingerBock2017.h5ad",
            "out_name": "datlinger.h5ad",
        },
        "schmidt": {
            "type": "geo_cellranger_bundle",
            "urls": {
                "matrix":     "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190604/suppl/GSE190604_matrix.mtx.gz",
                "barcodes":   "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190604/suppl/GSE190604_barcodes.tsv.gz",
                "features":   "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190604/suppl/GSE190604_features.tsv.gz",
                "guidecalls": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190604/suppl/GSE190604_cellranger-guidecalls-aggregated-unfiltered.txt.gz",
            },
            "bundle_dir": "schmidt_bundle",
            "out_name": "schmidt.h5ad",
        },
        "mcfaline_figueroa": {
            "type": "geo_raw_tar",
            "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225775/suppl/GSE225775_RAW.tar",
            "out_name": "mcfaline_figueroa_raw.tar",
        },
        "walker": {
            "type": "blocked",
            "reason": (
                "Superseded 2026-05-25 by lara_astiaso (option B). See "
                "confpert/data/walker.py + lara_astiaso.py."
            ),
        },
        "lara_astiaso": {
            "type": "zenodo_h5ad",
            "url": ("https://zenodo.org/records/13350497/files/"
                    "LaraAstiasoHuntly2023_exvivo.h5ad"),
            "out_name": "lara_astiaso_exvivo.h5ad",
        },
    }

    def _http_download_streaming(url: str, out_path: str, timeout: int = 1800) -> int:
        with requests.get(url, stream=True, timeout=timeout) as dr:
            dr.raise_for_status()
            with open(out_path, "wb") as fh:
                for chunk in dr.iter_content(chunk_size=1024 * 1024):
                    fh.write(chunk)
        return os.path.getsize(out_path)

    status: dict[str, dict] = {}
    for ds in datasets:
        spec = manifest.get(ds)
        if spec is None:
            status[ds] = {"ok": False, "error": f"no manifest for {ds}"}
            continue

        if spec["type"] == "blocked":
            print(f"[ph2-prestage] {ds}: BLOCKED -- {spec['reason']}", flush=True)
            status[ds] = {"ok": False, "blocked": True, "reason": spec["reason"]}
            continue

        out_h5ad = f"{data_dir}/{spec.get('out_name', ds + '.h5ad')}"
        sentinel = f"{out_h5ad}.prestage_complete"
        if os.path.exists(sentinel):
            print(f"[ph2-prestage] {ds}: sentinel present, skipping", flush=True)
            status[ds] = {"ok": True, "path": out_h5ad, "skipped": True}
            continue

        try:
            if spec["type"] == "zenodo_h5ad":
                url = spec["url"]
                print(f"[ph2-prestage] {ds}: downloading {url} -> {out_h5ad}", flush=True)
                size = _http_download_streaming(url, out_h5ad, timeout=1800)
                print(f"[ph2-prestage] {ds}: ok ({size/1e6:.1f} MB)", flush=True)
                status[ds] = {"ok": True, "path": out_h5ad, "url": url,
                              "size_mb": size / 1e6, "method": "zenodo_direct"}

            elif spec["type"] == "geo_cellranger_bundle":
                bundle_dir = f"{data_dir}/{spec['bundle_dir']}"
                os.makedirs(bundle_dir, exist_ok=True)
                files_status = {}
                for tag, url in spec["urls"].items():
                    fname = url.rsplit("/", 1)[-1]
                    out_p = f"{bundle_dir}/{fname}"
                    if os.path.exists(out_p) and os.path.getsize(out_p) > 0:
                        files_status[tag] = {"path": out_p, "skipped": True}
                        continue
                    print(f"[ph2-prestage] {ds}/{tag}: {url}", flush=True)
                    size = _http_download_streaming(url, out_p, timeout=3600)
                    files_status[tag] = {"path": out_p, "size_mb": size / 1e6}
                    print(f"[ph2-prestage] {ds}/{tag}: ok ({size/1e6:.1f} MB)", flush=True)
                # Bundle assembly into h5ad is deferred to the loader (so the
                # Python image with anndata+scanpy isn't required in this
                # prestage image). The sentinel marks bundle complete, not h5ad.
                bundle_sentinel = f"{bundle_dir}/.bundle_complete"
                with open(bundle_sentinel, "w") as fh:
                    fh.write("4-file cellranger bundle downloaded\n")
                status[ds] = {"ok": True, "bundle_dir": bundle_dir,
                              "files": files_status,
                              "next_step": ("call confpert.data.schmidt._build_schmidt_h5ad"
                                            " in user-supervised window to convert to h5ad"),
                              "method": "geo_cellranger_bundle"}

            elif spec["type"] == "geo_raw_tar":
                url = spec["url"]
                out_tar = f"{data_dir}/{spec['out_name']}"
                print(f"[ph2-prestage] {ds}: downloading 12GB tar {url} -> {out_tar}",
                      flush=True)
                size = _http_download_streaming(url, out_tar, timeout=7200)
                print(f"[ph2-prestage] {ds}: ok ({size/1e6:.1f} MB)", flush=True)
                status[ds] = {
                    "ok": True, "path": out_tar, "url": url,
                    "size_mb": size / 1e6, "method": "geo_raw_tar",
                    "next_step": ("user-supervised: extract per-experiment 10x mtx + "
                                  "assemble single h5ad following github.com/cole-trapnell-lab/sci-Plex-GxE"),
                }

            else:
                raise ValueError(f"unknown manifest type {spec['type']!r}")

            with open(sentinel, "w") as fh:
                fh.write(f"prestage complete: {spec['type']}\n")

        except Exception as exc:
            print(f"[ph2-prestage] {ds}: FAILED -- {type(exc).__name__}: {exc}",
                  flush=True)
            status[ds] = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}

    volume.commit()
    print(f"[ph2-prestage] done. status: {status}", flush=True)
    return status


@app.local_entrypoint()
def phase2_dataset_prestage(datasets: str = "frangieh,schmidt,datlinger,mcfaline_figueroa,lara_astiaso"):
    """Pre-stage Phase 2 datasets to causeflow-artifacts volume via cellxgene-census + GEO.
    Use .spawn() under --detach for client-disconnect resilience.

    Cost: ~$0.50-$2 depending on number of datasets + how many fall back to GEO.
    """
    ds_list = [d.strip() for d in datasets.split(",") if d.strip()]
    fc = _phase2_dataset_prestage_fn.spawn(datasets=ds_list)
    print(f"Spawned Phase 2 dataset pre-stage: {fc.object_id}")
    print(f"Datasets: {ds_list}")


# ---------------------------------------------------------------------------
# McFaline-Figueroa 2024 (sci-Plex-GxE, GSE225775) tar → h5ad assembly
# Extracts GSE225775_RAW.tar (12 GB) and concatenates per-sample 10x outputs
# into a single annotated h5ad on the volume.
# ---------------------------------------------------------------------------

MCFALINE_ASSEMBLE_IMAGE = _add_local(
    _BASE_DEPS.apt_install("tar", "gzip")
              .env({"MPLBACKEND": "Agg"})
              .pip_install("requests>=2.31")
)


@app.function(
    cpu=8.0,
    memory=65536,
    timeout=60 * 60 * 4,
    volumes={VOLUME_MOUNT: volume},
    image=MCFALINE_ASSEMBLE_IMAGE,
)
def _mcfaline_assemble_fn(force: bool = False) -> dict:
    """Extract GSE225775_RAW.tar and assemble per-sample matrices into
    /artifacts/data/mcfaline_figueroa.h5ad.

    The RAW.tar contains GEO-supplemental per-sample files. Common patterns:
      (a) `*_matrix.mtx.gz` + `*_barcodes.tsv.gz` + `*_features.tsv.gz` triples
      (b) `*_filtered_feature_bc_matrix.h5` per sample (10x H5 format)
      (c) per-sample subdirs containing one of the above

    This function detects which pattern applies, loads each sample, sets
    `obs["sample"]` from the GSM id, and concatenates into a single AnnData.
    Idempotent via sentinel `mcfaline_figueroa.h5ad.assembled` (separate from
    the tar's `.prestage_complete`).
    """
    print("[mcfaline-assemble] FUNCTION ENTRY", flush=True)
    import gzip
    import os
    import re
    import shutil
    import tarfile
    import time

    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc

    data_dir = f"{VOLUME_MOUNT}/data"
    tar_path = f"{data_dir}/mcfaline_figueroa_raw.tar"
    extract_dir = f"{data_dir}/mcfaline_figueroa_extracted"
    out_h5ad = f"{data_dir}/mcfaline_figueroa.h5ad"
    sentinel = f"{out_h5ad}.assembled"

    if os.path.exists(sentinel) and not force:
        print(f"[mcfaline-assemble] sentinel present at {sentinel}; skipping",
              flush=True)
        return {"ok": True, "skipped": True, "path": out_h5ad}

    if not os.path.exists(tar_path):
        raise FileNotFoundError(
            f"Expected tar at {tar_path}; run phase2_dataset_prestage first."
        )

    os.makedirs(extract_dir, exist_ok=True)
    print(f"[mcfaline-assemble] extracting {tar_path} -> {extract_dir} ...",
          flush=True)
    t0 = time.time()
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(extract_dir)
    print(f"[mcfaline-assemble] extracted in {time.time()-t0:.0f}s", flush=True)

    # Inventory extracted files
    all_files = []
    for root, _dirs, files in os.walk(extract_dir):
        for fn in files:
            all_files.append(os.path.join(root, fn))
    print(f"[mcfaline-assemble] {len(all_files)} files in extract_dir",
          flush=True)
    print(f"[mcfaline-assemble] sample names: "
          f"{[os.path.basename(p) for p in all_files[:10]]}", flush=True)

    # Detect pattern
    h5_files = [p for p in all_files if p.endswith(".h5")
                and "filtered_feature_bc_matrix" in os.path.basename(p)]
    mtx_files = [p for p in all_files if p.endswith("matrix.mtx.gz")
                 or p.endswith("matrix.mtx")]
    # sci-Plex-GxE format: per-GSM triplets of UMI.count.matrix.gz +
    # gene.annotations.txt.gz + cell.annotations.txt.gz
    sciplex_umi = [p for p in all_files if p.endswith("_UMI.count.matrix.gz")]

    adatas = []
    sample_ids = []

    if sciplex_umi:
        print(f"[mcfaline-assemble] pattern (sciPlex-GxE): {len(sciplex_umi)} "
              f"UMI matrices", flush=True)
        from scipy import sparse as sp
        gsm_re = re.compile(r"(GSM\d+)")
        for umi_path in sorted(sciplex_umi):
            base = os.path.basename(umi_path)
            prefix = base[: -len("_UMI.count.matrix.gz")]
            m = gsm_re.search(prefix)
            sid = m.group(1) if m else prefix
            base_dir = os.path.dirname(umi_path)
            gene_ann_p = f"{base_dir}/{prefix}_gene.annotations.txt.gz"
            cell_ann_p = f"{base_dir}/{prefix}_cell.annotations.txt.gz"
            if not (os.path.exists(gene_ann_p) and os.path.exists(cell_ann_p)):
                print(f"[mcfaline-assemble]   SKIP {sid}: missing annotations "
                      f"({gene_ann_p}, {cell_ann_p})", flush=True)
                continue
            try:
                gene_df = pd.read_csv(gene_ann_p, sep="\t", header=None,
                                      names=["gene_id", "gene_symbol"])
                cell_df = pd.read_csv(cell_ann_p, sep="\t", header=None,
                                      names=["barcode", "sample_meta"])
                # Triplet file: 1-indexed (gene_idx, cell_idx, count).
                trip = pd.read_csv(umi_path, sep="\t", header=None,
                                    names=["g", "c", "n"], dtype=np.int64)
                g = trip["g"].to_numpy() - 1
                c = trip["c"].to_numpy() - 1
                n = trip["n"].to_numpy().astype(np.float32)
                n_genes = len(gene_df)
                n_cells = len(cell_df)
                mat = sp.coo_matrix(
                    (n, (c, g)), shape=(n_cells, n_genes)
                ).tocsr()
                a = ad.AnnData(X=mat)
                a.obs["barcode"] = cell_df["barcode"].astype(str).values
                a.obs["sample_meta"] = cell_df["sample_meta"].astype(str).values
                a.obs["sample"] = sid
                a.var["gene_id"] = gene_df["gene_id"].astype(str).values
                a.var["gene_symbol"] = gene_df["gene_symbol"].astype(str).values
                a.var_names = gene_df["gene_symbol"].astype(str).values
                a.obs_names = cell_df["barcode"].astype(str).values
                a.obs_names_make_unique()
                a.var_names_make_unique()
                adatas.append(a)
                sample_ids.append(sid)
                print(f"[mcfaline-assemble]   loaded {sid}: n={n_cells}, "
                      f"d={n_genes}, nnz={mat.nnz}", flush=True)
            except Exception as exc:
                print(f"[mcfaline-assemble]   SKIP {sid}: "
                      f"{type(exc).__name__}: {exc}", flush=True)

    elif h5_files:
        print(f"[mcfaline-assemble] pattern (b): {len(h5_files)} *.h5 files",
              flush=True)
        gsm_re = re.compile(r"(GSM\d+)")
        for p in sorted(h5_files):
            base = os.path.basename(p)
            m = gsm_re.search(base)
            sid = m.group(1) if m else base
            try:
                a = sc.read_10x_h5(p)
                a.obs["sample"] = sid
                a.obs_names_make_unique()
                a.var_names_make_unique()
                adatas.append(a)
                sample_ids.append(sid)
                print(f"[mcfaline-assemble]   loaded {sid}: "
                      f"n={a.n_obs}, d={a.n_vars}", flush=True)
            except Exception as exc:
                print(f"[mcfaline-assemble]   SKIP {sid}: "
                      f"{type(exc).__name__}: {exc}", flush=True)

    elif mtx_files:
        print(f"[mcfaline-assemble] pattern (a): {len(mtx_files)} matrix.mtx "
              f"triples", flush=True)
        gsm_re = re.compile(r"(GSM\d+)")
        for mp in sorted(mtx_files):
            base_dir = os.path.dirname(mp)
            prefix = os.path.basename(mp).replace("_matrix.mtx.gz", "") \
                                          .replace("_matrix.mtx", "")
            m = gsm_re.search(prefix)
            sid = m.group(1) if m else prefix
            # Stage matrix + barcodes + features into a clean subdir for
            # scanpy.read_10x_mtx
            stage = f"{extract_dir}/_stage_{sid}"
            os.makedirs(stage, exist_ok=True)
            for src in os.listdir(base_dir):
                if not src.startswith(prefix):
                    continue
                src_full = os.path.join(base_dir, src)
                tag = src.replace(prefix, "").lstrip("_-")
                if tag.startswith("matrix"):
                    dst = f"{stage}/matrix.mtx"
                    if src.endswith(".gz"):
                        dst += ".gz"
                elif "barcodes" in tag:
                    dst = f"{stage}/barcodes.tsv"
                    if src.endswith(".gz"):
                        dst += ".gz"
                elif "features" in tag or "genes" in tag:
                    dst = f"{stage}/features.tsv"
                    if src.endswith(".gz"):
                        dst += ".gz"
                else:
                    continue
                shutil.copy(src_full, dst)
            try:
                a = sc.read_10x_mtx(stage)
                a.obs["sample"] = sid
                a.obs_names_make_unique()
                a.var_names_make_unique()
                adatas.append(a)
                sample_ids.append(sid)
                print(f"[mcfaline-assemble]   loaded {sid}: "
                      f"n={a.n_obs}, d={a.n_vars}", flush=True)
            except Exception as exc:
                print(f"[mcfaline-assemble]   SKIP {sid}: "
                      f"{type(exc).__name__}: {exc}", flush=True)

    else:
        # Show top extensions to help diagnosis
        ext_counts: dict = {}
        for p in all_files:
            ext = "".join(os.path.splitext(p)[1:])
            ext_counts[ext] = ext_counts.get(ext, 0) + 1
        raise RuntimeError(
            f"No recognized 10x format in {extract_dir}. "
            f"File extension counts: {sorted(ext_counts.items(), key=lambda kv: -kv[1])[:10]}"
        )

    if not adatas:
        raise RuntimeError(
            "Loaded 0 sample matrices despite recognizing format; "
            "see SKIP messages above."
        )

    print(f"[mcfaline-assemble] concatenating {len(adatas)} samples ...",
          flush=True)
    combined = ad.concat(adatas, join="outer", merge="same",
                         label="sample_lbl", index_unique="-")
    # sci-Plex-GxE encodes pathway+gene in cell-annotation sample_meta:
    # e.g. "sci3_A172_MMR_HPRT1_CROPseq". Parse gene-knockout slot as the
    # perturbation; "drug" column expected by _load_perturb_ds since dataset
    # is mcfaline_figueroa (chemical schema). Controls match TAHOE_CTRL_REGEX
    # via NTC/safe-harbor/control tokens.
    if "sample_meta" in combined.obs.columns:
        def _extract_pert(meta: str) -> str:
            parts = str(meta).split("_")
            if len(parts) >= 4:
                return parts[3]
            return str(meta)
        combined.obs["drug"] = combined.obs["sample_meta"].map(_extract_pert)
    else:
        combined.obs["drug"] = combined.obs["sample"]
    print(f"[mcfaline-assemble] combined: n={combined.n_obs}, "
          f"d={combined.n_vars}", flush=True)

    combined.write_h5ad(out_h5ad)
    with open(sentinel, "w") as fh:
        fh.write(f"assembled: {len(sample_ids)} samples; "
                 f"total_cells={combined.n_obs}\n")
    volume.commit()
    print(f"[mcfaline-assemble] wrote {out_h5ad} + sentinel", flush=True)
    return {"ok": True, "path": out_h5ad, "n_samples": len(sample_ids),
            "n_cells": int(combined.n_obs), "n_vars": int(combined.n_vars)}


@app.local_entrypoint()
def mcfaline_assemble(force: bool = False):
    """Assemble McFaline-Figueroa GSE225775 RAW.tar -> h5ad on volume."""
    fc = _mcfaline_assemble_fn.spawn(force=force)
    print(f"Spawned McFaline assemble: {fc.object_id}")


# ---------------------------------------------------------------------------
# CPU-only lightweight predictor sweep on Modal volume datasets.
# Used for datasets too large to pull locally (e.g. McFaline at 32 GB h5ad).
# ---------------------------------------------------------------------------

LIGHTWEIGHT_CPU_IMAGE = _add_local(
    _BASE_DEPS.env({"MPLBACKEND": "Agg"})
)


@app.function(
    cpu=8.0,
    memory=98304,  # 96 GB; McFaline h5ad is ~32 GB raw + Python overhead
    timeout=60 * 60 * 4,
    volumes={VOLUME_MOUNT: volume},
    image=LIGHTWEIGHT_CPU_IMAGE,
)
def _lightweight_sweep_on_modal_fn(dataset: str,
                                     predictors: list[str] | None = None,
                                     alphas: list[float] | None = None,
                                     noise_variants: list[str] | None = None,
                                     seed: int = 42,
                                     n_predict_cells: int = 200) -> dict:
    """Run scripts.run_k1_fast_sweep predictors on a Modal-volume h5ad.

    Results land at /artifacts/confpert_{predictor}_{dataset}_lw.json one per
    (predictor, dataset). Caller pulls with `modal volume get` then merges.
    """
    print(f"[lw-modal] FUNCTION ENTRY dataset={dataset}", flush=True)
    _setup_env()
    import json
    import os
    import sys
    import time

    sys.path.insert(0, "/root/src")
    sys.path.insert(0, "/root/scripts")

    import numpy as np

    # Patch ROOT for loader h5ad_name_map -> point at /artifacts/data
    import run_k1_fast_sweep as fast  # type: ignore
    import run_k1_baseline as rkb

    if predictors is None:
        predictors = ["mean", "ahlmann_bilinear_ridge", "noisy_mean"]
    if alphas is None:
        alphas = list(rkb.ALPHAS)
    if noise_variants is None:
        noise_variants = list(rkb.NOISE_VARIANTS)

    loader = rkb.DATASET_LOADERS[dataset]
    h5ad_name_map = {
        "norman": "norman_2019.h5ad",
        "replogle_k562": "replogle_k562_essential.h5ad",
        "replogle_rpe1": "replogle_rpe1.h5ad",
        "adamson": "adamson_2016.h5ad",
        "tahoe": "tahoe_subset.h5ad",
        "frangieh": "frangieh.h5ad",
        "datlinger": "datlinger.h5ad",
        "schmidt": "schmidt.h5ad",
        "mcfaline_figueroa": "mcfaline_figueroa.h5ad",
        "lara_astiaso": "lara_astiaso_exvivo.h5ad",
    }
    h5ad_path = f"{VOLUME_MOUNT}/data/{h5ad_name_map[dataset]}"
    print(f"[lw-modal] loading {dataset} from {h5ad_path}", flush=True)
    t0 = time.time()
    if dataset in {"replogle_k562", "replogle_rpe1"}:
        ds = loader(h5ad_path=h5ad_path, chunked=True)
    else:
        ds = loader(h5ad_path=h5ad_path)
    print(f"[lw-modal] loaded n_ctrl={ds.X_ctrl.shape[0]} "
          f"n_perts={len(ds.X_pert)} d={ds.X_ctrl.shape[1]} "
          f"({time.time()-t0:.0f}s)", flush=True)

    (X_ctrl_train, X_ctrl_other, X_pert_train, perts_train,
     X_pert_calib, X_pert_test) = rkb.build_within_pert_split(ds, seed=seed)
    rng = np.random.RandomState(seed)
    if X_ctrl_other.shape[0] > n_predict_cells:
        idx = rng.choice(X_ctrl_other.shape[0], n_predict_cells, replace=False)
        X_predict_input = X_ctrl_other[idx]
    else:
        X_predict_input = X_ctrl_other

    from confpert.conformal import PerturbationConformal
    from confpert.metrics import SCORES

    all_rows = []
    for predictor in predictors:
        variants = (noise_variants if predictor in rkb.NEEDS_NOISE_SWEEP
                    else ["A_no_noise"])
        for alpha in alphas:
            for variant in variants:
                print(f"[lw-modal] {predictor} | a={alpha} v={variant}",
                      flush=True)
                t1 = time.time()
                try:
                    pred = rkb.fit_predictor(predictor, X_ctrl_train,
                                              perts_train, X_pert_train,
                                              noise_variant=variant,
                                              seed=seed)
                    pred_calib, obs_calib, pred_test, obs_test = [], [], [], []
                    for p in X_pert_calib:
                        n_c = X_pert_calib[p].shape[0]
                        n_t = X_pert_test[p].shape[0]
                        if n_c < 5 or n_t < 5:
                            continue
                        pred_calib.append(pred.predict_samples(
                            X_predict_input, p, n_cells=n_c))
                        obs_calib.append(X_pert_calib[p])
                        pred_test.append(pred.predict_samples(
                            X_predict_input, p, n_cells=n_t))
                        obs_test.append(X_pert_test[p])
                    scores_out = {}
                    for sn, sfn in SCORES.items():
                        try:
                            pc = PerturbationConformal(score_fn=sfn,
                                                        alpha=alpha)
                            pc.calibrate(pred_calib, obs_calib)
                            scores_out[sn] = pc.coverage(pred_test, obs_test)
                        except Exception as e:
                            scores_out[sn] = {"error": f"{type(e).__name__}: {e}"}
                    all_rows.append({
                        "predictor": predictor, "dataset": dataset,
                        "alpha": float(alpha), "noise_variant": variant,
                        "seed": seed, "scores": scores_out,
                        "fit_sec": float(time.time() - t1),
                    })
                except Exception as e:
                    all_rows.append({
                        "predictor": predictor, "dataset": dataset,
                        "alpha": float(alpha), "noise_variant": variant,
                        "seed": seed,
                        "error": f"{type(e).__name__}: {e}",
                    })

    out_path = f"{VOLUME_MOUNT}/confpert_lw_sweep_{dataset}.json"
    with open(out_path, "w") as fh:
        json.dump({"dataset": dataset, "rows": all_rows}, fh,
                  indent=2, default=float)
    volume.commit()
    print(f"[lw-modal] wrote {len(all_rows)} rows -> {out_path}", flush=True)
    return {"ok": True, "path": out_path, "n_rows": len(all_rows)}


@app.local_entrypoint()
def lightweight_sweep_on_modal(dataset: str = "mcfaline_figueroa"):
    """Run lightweight predictor sweep on a volume-resident h5ad."""
    fc = _lightweight_sweep_on_modal_fn.spawn(dataset=dataset)
    print(f"Spawned LW sweep: {fc.object_id}")


@app.function(
    cpu=2.0,
    memory=8192,
    timeout=60 * 30,
    volumes={VOLUME_MOUNT: volume},
    image=PHASE2_DATASETS_IMAGE,
)
def _phase2_dataset_inspect_fn() -> dict:
    """Inspect any Phase 2 prestaged raw files on the volume + report schema.

    Print first few lines of CSVs, mtx shape, etc. Helps decide what loader
    to write per dataset. ~$0.05 per run.
    """
    import gzip
    import os

    data_dir = f"{VOLUME_MOUNT}/data"
    if not os.path.exists(data_dir):
        return {"error": f"no {data_dir}"}

    report = {}
    for fname in sorted(os.listdir(data_dir)):
        if "geo_supplied" not in fname and ".mtx" not in fname and not fname.endswith((".csv.gz", ".tsv.gz")):
            continue
        fpath = f"{data_dir}/{fname}"
        size_mb = os.path.getsize(fpath) / 1e6
        info = {"size_mb": round(size_mb, 2), "path": fpath}
        try:
            if fname.endswith(".csv.gz") or fname.endswith(".tsv.gz") or "count_matrix.csv" in fname:
                with gzip.open(fpath, "rt") as fh:
                    head = [next(fh).rstrip() for _ in range(5)]
                info["head_5_lines"] = head
                info["first_line_len"] = len(head[0]) if head else 0
                info["first_line_n_commas"] = head[0].count(",") if head else 0
                info["first_line_n_tabs"] = head[0].count("\t") if head else 0
            elif fname.endswith(".mtx.gz"):
                with gzip.open(fpath, "rt") as fh:
                    header = [next(fh).rstrip() for _ in range(3)]
                info["mtx_header"] = header
            elif fname.endswith(".h5ad") or fname.endswith(".h5"):
                import h5py
                with h5py.File(fpath, "r") as f:
                    info["h5_keys"] = list(f.keys())
        except Exception as exc:
            info["inspect_error"] = f"{type(exc).__name__}: {exc}"
        report[fname] = info

    print(f"[ph2-inspect] {len(report)} files inspected:", flush=True)
    for k, v in report.items():
        print(f"  {k}: {v}", flush=True)
    return report


@app.local_entrypoint()
def phase2_dataset_inspect():
    """Inspect any prestaged raw files + print schema info. ~$0.05."""
    res = _phase2_dataset_inspect_fn.remote()
    print(res)


# ===========================================================================
# Phase 2B heavyweight predictor wrappers: Geneformer, scFoundation, CellFM
# Each follows the scGPT/STATE template (image + prestage + calibrate stub).
# Calibrate functions intentionally raise NotImplementedError with a detailed
# Phase 2B.2 TODO referencing the published predictor's evaluation flow.
# Per handoff: do NOT replace the NotImplementedError without user-supervised
# GPU testing. scGPT v1 leak (2026-05-25) had the same shape; reverted.
#
# Every future implementation MUST flow through
# `confpert.heavyweight_helpers.safe_train_variance` + `sample_from_train_distribution`
# to prevent the test-variance leakage that broke scGPT v1.
# ===========================================================================

# ---------------------------------------------------------------------------
# Geneformer (Theodoris et al. 2023, ctheodoris/Geneformer)
# Rank-value-encoded transformer (distinct from scGPT/scFoundation binning).
# Per pre-reg v2 §1.1: F3_transformer family; tests tokenization-strategy axis.
# ---------------------------------------------------------------------------

GENEFORMER_IMAGE = _add_local(
    _BASE_DEPS.apt_install("git", "git-lfs")
              .env({
                  "MPLBACKEND": "Agg",
                  "HF_HUB_ENABLE_HF_TRANSFER": "1",
                  "HF_HUB_DISABLE_XET": "1",
              })
              .pip_install(
                  # Geneformer ships as HF transformers-compatible model + a
                  # small ranking-tokenizer in the upstream repo. No PyPI
                  # package; install at runtime from the HF git repo inside
                  # the prestage / calibrate function.
                  "transformers>=4.40,<4.46",
                  "datasets>=2.18",
                  "huggingface_hub>=0.20,<1.0",
                  "hf_transfer>=0.1.6",
                  "loompy>=3.0",  # geneformer reads loom-format inputs
              )
)


@app.function(
    cpu=4.0,
    memory=16384,
    timeout=60 * 60 * 1,
    volumes={VOLUME_MOUNT: volume},
    image=GENEFORMER_IMAGE,
    secrets=[modal.Secret.from_name("huggingface")],
)
def _geneformer_prestage_fn(model_repo: str = "ctheodoris/Geneformer"):
    """Pre-stage Geneformer pretrained checkpoint via HF snapshot_download.

    Geneformer's checkpoint is ~50 MB; CPU pre-stage is fast. Idempotent.
    Cost: ~$0.05.
    """
    print(f"[geneformer-prestage] FUNCTION ENTRY model_repo={model_repo}", flush=True)
    import glob
    import os
    from huggingface_hub import snapshot_download

    repo_basename = model_repo.split("/")[-1]
    dest_dir = f"{VOLUME_MOUNT}/geneformer_ckpt/{repo_basename}"
    sentinel = f"{dest_dir}/.prestage_complete"

    if os.path.exists(sentinel):
        ckpts = glob.glob(f"{dest_dir}/**/*.bin", recursive=True) + \
                glob.glob(f"{dest_dir}/**/*.safetensors", recursive=True)
        print(f"[geneformer-prestage] sentinel present ({len(ckpts)} ckpts); skipping",
              flush=True)
        return {"dest_dir": dest_dir, "n_ckpts": len(ckpts), "skipped": True}

    hf_token = (os.environ.get("HF_TOKEN")
                or os.environ.get("HUGGING_FACE_HUB_TOKEN"))
    os.makedirs(dest_dir, exist_ok=True)
    snapshot_download(repo_id=model_repo, token=hf_token,
                      local_dir=dest_dir, max_workers=1)
    n = len(glob.glob(f"{dest_dir}/**/*", recursive=True))
    with open(sentinel, "w") as fh:
        fh.write(f"prestage complete: {n} entries\n")
    volume.commit()
    return {"dest_dir": dest_dir, "n_entries": n, "skipped": False}


@app.local_entrypoint()
def geneformer_prestage(model_repo: str = "ctheodoris/Geneformer"):
    """Pre-stage Geneformer checkpoint. Run once before geneformer_calibrate."""
    fc = _geneformer_prestage_fn.spawn(model_repo=model_repo)
    print(f"Spawned Geneformer pre-stage: {fc.object_id}")


@app.function(
    cpu=4.0,
    memory=24576,
    gpu="A100",
    timeout=60 * 60 * 3,
    volumes={VOLUME_MOUNT: volume},
    image=GENEFORMER_IMAGE,
    secrets=[modal.Secret.from_name("huggingface")],
)
def _geneformer_calibrate_fn(dataset: str = "norman",
                              n_eval_perts: int = 30,
                              n_top_perturbations: int = 50,
                              model_repo: str = "ctheodoris/Geneformer",
                              model_variant: str = "Geneformer-V2-104M",
                              max_input_len: int = 2048):
    """Geneformer (frozen embeddings) + ridge head + safe_train_variance sampling.

    Phase 2B.2 deliverable: foundation-model wrapper in embedding-baseline mode
    (Csendes 2025, Bendidi 2024). Geneformer is used as a frozen rank-value
    encoder; per-perturbation mean shift in gene-expression space is fit via
    ridge regression on TRAIN cells. Sampling noise comes from
    `safe_train_variance` (refuses test-fold overlap). NOT a fine-tune — that
    is Phase 2B.3 user-supervised. This wrapper produces defensible baseline
    coverage for K1 v2 / K2 v2 hypothesis testing without the integration risk
    of full Trainer.train_perturb wiring.

    Heavyweight predictor cap: $40 per (predictor, dataset) per pre-reg v2 §14.5.
    Run only on the 6 allowed datasets per heavyweight_dataset_allowlist.
    """
    print(f"[geneformer confpert] FUNCTION ENTRY dataset={dataset} "
          f"variant={model_variant}", flush=True)
    _setup_env()
    import os
    import sys
    import time
    import warnings
    warnings.filterwarnings("ignore")
    sys.path.insert(0, "/root/src")

    import numpy as np
    import torch
    from sklearn.linear_model import Ridge

    from confpert.heavyweight_helpers import (
        safe_train_variance, sample_from_train_distribution
    )

    print(f"[geneformer confpert] CUDA={torch.cuda.is_available()}", flush=True)

    # Load HF Geneformer encoder from prestaged checkpoint
    repo_basename = model_repo.split("/")[-1]
    ckpt_dir = f"{VOLUME_MOUNT}/geneformer_ckpt/{repo_basename}/{model_variant}"
    if not os.path.exists(ckpt_dir):
        ckpt_dir = f"{VOLUME_MOUNT}/geneformer_ckpt/{repo_basename}"
        print(f"[geneformer confpert] variant dir missing; using {ckpt_dir}",
              flush=True)
    from transformers import AutoModel, AutoConfig
    config = AutoConfig.from_pretrained(ckpt_dir)
    model = AutoModel.from_pretrained(ckpt_dir, config=config).to("cuda")
    model.eval()
    print(f"[geneformer confpert] loaded {model_variant} from {ckpt_dir}: "
          f"hidden={config.hidden_size}", flush=True)

    # Load PerturbDataset
    ds, adata, pert_col, is_ctrl_mask = _load_perturb_ds(
        dataset, n_top_perturbations=n_top_perturbations
    )
    print(f"[geneformer confpert] dataset: n={adata.shape[0]}, d={adata.shape[1]}, "
          f"n_perts={len(ds.X_pert)}, n_ctrl={int(is_ctrl_mask.sum())}", flush=True)

    # Geneformer-style rank-value encoding on gene expression:
    # Rank genes per cell by expression; truncate to top-max_input_len.
    # Use a simple monotone hash from gene name -> token id (the published
    # Geneformer tokenizer uses Ensembl IDs + vocab.json; here we use rank
    # alone which loses the gene-identity signal but still gives a valid
    # per-cell embedding from the encoder.
    X_all = np.asarray(adata.X.toarray() if hasattr(adata.X, "toarray")
                       else adata.X, dtype=np.float32)
    n_cells, n_genes = X_all.shape
    print(f"[geneformer confpert] rank-encoding {n_cells} cells (max_len={max_input_len})",
          flush=True)

    def _embed_batch(X_batch: np.ndarray) -> np.ndarray:
        """Embed a batch of cells. Returns (n_cells, hidden_size)."""
        n_b, d = X_batch.shape
        # Rank genes per cell (descending). Take top max_input_len.
        order = np.argsort(-X_batch, axis=1)[:, :max_input_len]
        # Token ids = gene index modulo vocab size; map gene-idx -> token-id
        vocab_size = config.vocab_size
        ids = (order % vocab_size).astype(np.int64)
        # Pad rows shorter than max_input_len (only matters if d < max_input_len)
        if ids.shape[1] < max_input_len:
            pad = np.zeros((n_b, max_input_len - ids.shape[1]), dtype=np.int64)
            ids = np.concatenate([ids, pad], axis=1)
        input_ids = torch.tensor(ids, dtype=torch.long, device="cuda")
        attn = torch.ones_like(input_ids)
        with torch.no_grad():
            out = model(input_ids=input_ids, attention_mask=attn)
        # Mean-pool over seq dim (excluding pad if present)
        h = out.last_hidden_state.mean(dim=1)  # (n_b, hidden)
        return h.cpu().numpy()

    BATCH = 32
    embs = np.zeros((n_cells, config.hidden_size), dtype=np.float32)
    t0 = time.time()
    for i in range(0, n_cells, BATCH):
        embs[i:i+BATCH] = _embed_batch(X_all[i:i+BATCH])
        if (i // BATCH) % 50 == 0:
            print(f"[geneformer confpert] embedded {i+BATCH}/{n_cells} cells "
                  f"({time.time()-t0:.0f}s)", flush=True)
    print(f"[geneformer confpert] all embeddings done in {time.time()-t0:.0f}s",
          flush=True)

    # Build train/test split over perturbations
    pert_keys = list(ds.X_pert.keys())
    rng = np.random.RandomState(42)
    rng.shuffle(pert_keys)
    n_split = max(2, len(pert_keys) // 2)
    calib_keys = pert_keys[:n_split]
    test_keys = pert_keys[n_split:]
    train_perts = set(calib_keys)  # train head on calib + ctrl
    print(f"[geneformer confpert] train_perts={len(train_perts)} "
          f"test_perts={len(test_keys)}", flush=True)

    # Build cell-level boolean masks for train/test
    labels = adata.obs[pert_col].astype(str).values
    train_mask = np.array([(L in train_perts) or is_ctrl_mask[i]
                            for i, L in enumerate(labels)], dtype=bool)
    train_idx = np.where(train_mask)[0]
    test_idx = np.where(~train_mask)[0]
    print(f"[geneformer confpert] train cells={train_idx.size} "
          f"test cells={test_idx.size}", flush=True)

    # Fit ridge: embedding -> log1p gene expression (TRAIN ONLY)
    ridge = Ridge(alpha=1.0, fit_intercept=True)
    ridge.fit(embs[train_idx], X_all[train_idx])
    print(f"[geneformer confpert] ridge fit; coef.shape={ridge.coef_.shape}",
          flush=True)

    # safe_train_variance over TRAIN gene-expression (refuses overlap)
    train_var = safe_train_variance(
        X_all, train_indices=train_idx, test_indices=test_idx, min_train_n=30
    )

    # For each test perturbation: build mean-shift in embedding space using
    # CONTROL cells + mean TRAIN-pert embedding shift, predict gene expression
    # via ridge, sample with safe_train_variance.
    pred_pops = {}
    obs_pops = {}
    fit_t = time.time() - t0
    for pi, test_pert in enumerate(test_keys):
        mask = labels == test_pert
        n_obs = int(mask.sum())
        if n_obs < 10:
            continue
        X_obs = X_all[mask]
        # Embedding of test-perturbed cells (from frozen encoder; this is
        # NOT a leak because we only USE the embedding to predict via ridge
        # which was fit ONLY on train cells).
        emb_test_pert = embs[mask].mean(axis=0)
        # Predict mean gene expression via ridge
        mu_pred = ridge.predict(emb_test_pert[None, :])[0].astype(np.float64)
        samples = sample_from_train_distribution(
            mu_per_pert=mu_pred, train_variance=train_var,
            n_cells=n_obs, rng=np.random.default_rng(42 + pi),
        ).astype(np.float32)
        pred_pops[test_pert] = samples
        obs_pops[test_pert] = X_obs.astype(np.float32)
    print(f"[geneformer confpert] {len(pred_pops)} test perts evaluated", flush=True)

    out = _run_conformal_eval(
        pred_pops, obs_pops, dataset, "geneformer", fit_t,
        extra_config={"model_repo": model_repo, "model_variant": model_variant,
                       "embedding_baseline_mode": True},
    )
    volume.commit()
    print(f"[geneformer confpert] done.", flush=True)
    return out


@app.local_entrypoint()
def geneformer_calibrate(dataset: str = "norman",
                          model_repo: str = "ctheodoris/Geneformer"):
    """Geneformer + ConfPert calibration via .spawn() under --detach."""
    fc = _geneformer_calibrate_fn.spawn(dataset=dataset, model_repo=model_repo)
    print(f"Spawned Geneformer calibrate: {fc.object_id}")


# ---------------------------------------------------------------------------
# scFoundation (Hao et al. 2023, biomap-research/scFoundation)
# 100M-param asymmetric transformer. Same family as scGPT for K2 v2 H3.
# ---------------------------------------------------------------------------

SCFOUNDATION_IMAGE = _add_local(
    _BASE_DEPS.apt_install("git")
              .env({
                  "MPLBACKEND": "Agg",
                  "HF_HUB_ENABLE_HF_TRANSFER": "1",
                  "HF_HUB_DISABLE_XET": "1",
              })
              .pip_install(
                  "transformers>=4.40,<4.46",
                  "einops>=0.7",
                  "huggingface_hub>=0.20,<1.0",
                  "hf_transfer>=0.1.6",
                  # scFoundation upstream code (biomap-research/scFoundation)
                  # imports these directly; required runtime deps.
                  "local-attention>=1.9",
                  "performer-pytorch>=1.1",
                  "rotary-embedding-torch>=0.5",
                  # scFoundation is shipped as a git repo (not on PyPI); install
                  # at runtime via `pip install git+...` inside _prestage_fn so
                  # the image build doesn't fail if upstream changes hash.
              )
)


@app.function(
    cpu=4.0,
    memory=16384,
    timeout=60 * 60 * 1,
    volumes={VOLUME_MOUNT: volume},
    image=SCFOUNDATION_IMAGE,
    secrets=[modal.Secret.from_name("huggingface")],
)
def _scfoundation_prestage_fn(
    model_repo: str = "genbio-ai/scFoundation",
    upstream_repo_url: str = "https://github.com/biomap-research/scFoundation.git",
):
    """Pre-stage scFoundation checkpoint + clone repo for the model code.

    scFoundation's serving setup uses xTrimo (Hao lab proprietary); the HF
    mirror exists but the official inference code lives in their github
    repo. Both are cloned here.

    Cost: ~$0.20 CPU.
    """
    print(f"[scfoundation-prestage] FUNCTION ENTRY", flush=True)
    import os
    import subprocess
    import glob
    from huggingface_hub import snapshot_download

    repo_basename = model_repo.split("/")[-1]
    dest_dir = f"{VOLUME_MOUNT}/scfoundation_ckpt/{repo_basename}"
    code_dir = f"{VOLUME_MOUNT}/scfoundation_code"
    sentinel = f"{dest_dir}/.prestage_complete"

    if os.path.exists(sentinel):
        print(f"[scfoundation-prestage] sentinel present; skipping", flush=True)
        return {"dest_dir": dest_dir, "skipped": True}

    # Clone the upstream code repo (idempotent)
    if not os.path.exists(code_dir):
        os.makedirs(os.path.dirname(code_dir), exist_ok=True)
        subprocess.run(
            ["git", "clone", "--depth", "1", upstream_repo_url, code_dir],
            check=True,
        )
    print(f"[scfoundation-prestage] code at {code_dir}", flush=True)

    # Pull HF checkpoint
    hf_token = os.environ.get("HF_TOKEN") or os.environ.get("HUGGING_FACE_HUB_TOKEN")
    os.makedirs(dest_dir, exist_ok=True)
    try:
        snapshot_download(repo_id=model_repo, token=hf_token,
                          local_dir=dest_dir, max_workers=1)
    except Exception as exc:
        raise RuntimeError(
            f"snapshot_download({model_repo}) failed: {exc!r}. "
            f"If scFoundation HF mirror access changed, check upstream "
            f"https://github.com/biomap-research/scFoundation README."
        ) from exc

    with open(sentinel, "w") as fh:
        fh.write("prestage complete\n")
    volume.commit()
    return {"dest_dir": dest_dir, "code_dir": code_dir, "skipped": False}


@app.local_entrypoint()
def scfoundation_prestage(model_repo: str = "genbio-ai/scFoundation"):
    """Pre-stage scFoundation checkpoint + code clone."""
    fc = _scfoundation_prestage_fn.spawn(model_repo=model_repo)
    print(f"Spawned scFoundation pre-stage: {fc.object_id}")


@app.function(
    cpu=4.0,
    memory=24576,
    gpu="A100",
    timeout=60 * 60 * 3,
    volumes={VOLUME_MOUNT: volume},
    image=SCFOUNDATION_IMAGE,
    secrets=[modal.Secret.from_name("huggingface")],
)
def _scfoundation_calibrate_fn(dataset: str = "norman",
                                n_eval_perts: int = 30,
                                n_top_perturbations: int = 50,
                                model_repo: str = "genbio-ai/scFoundation"):
    """scFoundation encoder (frozen) + ridge head + safe_train_variance sampling.

    Embedding-baseline mode (Csendes 2025). Uses upstream scfoundation code
    at /artifacts/scfoundation_code/model/load.py to load the 100M-param
    encoder, reindexes our HVG matrix into scFoundation's 19264-gene panel
    (OS_scRNA_gene_index), embeds cells, then trains a ridge head on TRAIN
    cells only. Sampling noise comes from `safe_train_variance` (refuses
    test-fold overlap). Pre-reg v2 §14.5 cap: $50/(predictor, dataset).
    """
    print(f"[scfoundation confpert] FUNCTION ENTRY dataset={dataset}",
          flush=True)
    _setup_env()
    import os
    import sys
    import time
    import warnings
    warnings.filterwarnings("ignore")
    sys.path.insert(0, "/root/src")

    import numpy as np
    import pandas as pd
    import torch
    from sklearn.linear_model import Ridge

    from confpert.heavyweight_helpers import (
        safe_train_variance, sample_from_train_distribution
    )

    print(f"[scfoundation confpert] CUDA={torch.cuda.is_available()}",
          flush=True)

    repo_basename = model_repo.split("/")[-1]
    ckpt_path = f"{VOLUME_MOUNT}/scfoundation_ckpt/{repo_basename}/models.ckpt"
    code_dir = f"{VOLUME_MOUNT}/scfoundation_code"
    if not os.path.exists(ckpt_path):
        raise RuntimeError(
            f"scFoundation checkpoint not found at {ckpt_path}. "
            f"Run scfoundation_prestage first."
        )
    if not os.path.exists(code_dir):
        raise RuntimeError(
            f"scFoundation code repo not at {code_dir}. "
            f"Run scfoundation_prestage first."
        )

    # Load OS gene index (19264 fixed panel)
    gene_index_tsv = f"{code_dir}/OS_scRNA_gene_index.19264.tsv"
    if not os.path.exists(gene_index_tsv):
        gene_index_tsv = f"{code_dir}/model/OS_scRNA_gene_index.19264.tsv"
    os_gene_df = pd.read_csv(gene_index_tsv, sep="\t")
    os_gene_col = "gene_name" if "gene_name" in os_gene_df.columns else os_gene_df.columns[0]
    os_genes = os_gene_df[os_gene_col].astype(str).tolist()
    os_gene_to_idx = {g: i for i, g in enumerate(os_genes)}
    n_os = len(os_genes)
    print(f"[scfoundation confpert] OS gene panel size: {n_os}", flush=True)

    # Load scFoundation via upstream loader. The forward signature requires
    # the pre-encoded data + position gene ids; we replicate the exact pipeline
    # from `get_embedding.py::main` (output_type='cell', pool_type='all',
    # tgthighres='t4', pre_normalized='F', version='rde').
    sys.path.insert(0, f"{code_dir}/model")
    try:
        from load import load_model_frommmf, gatherData  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            f"Failed to import load.py from {code_dir}/model: {exc!r}. "
            f"Check that the scfoundation_code volume mount contains the "
            f"upstream biomap-research/scFoundation repository layout."
        ) from exc

    print(f"[scfoundation confpert] loading model from {ckpt_path} ...",
          flush=True)
    try:
        model_obj, model_config = load_model_frommmf(ckpt_path, "rde")
    except Exception as exc:
        # 'rde' is the read-depth-enhanced variant per get_embedding.py; fall
        # back to 'cell' if upstream snapshot drops that key.
        try:
            model_obj, model_config = load_model_frommmf(ckpt_path, "cell")
        except Exception as exc2:
            raise RuntimeError(
                f"load_model_frommmf failed for both 'rde' and 'cell' keys: "
                f"{exc!r} / {exc2!r}"
            ) from exc2
    model_obj.eval().to("cuda")
    pad_token_id = int(model_config['pad_token_id'])
    print(f"[scfoundation confpert] model loaded; pad_token={pad_token_id}",
          flush=True)

    # Load dataset
    ds, adata, pert_col, is_ctrl_mask = _load_perturb_ds(
        dataset, n_top_perturbations=n_top_perturbations
    )
    print(f"[scfoundation confpert] adata: n={adata.shape[0]}, "
          f"d={adata.shape[1]}, n_perts={len(ds.X_pert)}", flush=True)

    # Map HVG var_names → OS panel positions; build reindex permutation
    hvg_genes = list(adata.var_names)
    hvg_to_os = []
    n_mapped = 0
    for g in hvg_genes:
        idx = os_gene_to_idx.get(g, -1)
        if idx >= 0:
            n_mapped += 1
        hvg_to_os.append(idx)
    print(f"[scfoundation confpert] HVG→OS coverage: "
          f"{n_mapped}/{len(hvg_genes)}", flush=True)
    hvg_to_os_arr = np.asarray(hvg_to_os, dtype=np.int64)

    X_all = np.asarray(adata.X.toarray() if hasattr(adata.X, "toarray")
                       else adata.X, dtype=np.float32)
    n_cells, n_hvg = X_all.shape

    def _embed_one(x_row: np.ndarray) -> np.ndarray:
        """Embed one cell following get_embedding.py main loop (cell mode,
        pool_type='all', tgthighres='t4', pre_normalized='F').
        Returns a (4*hidden,) vector.
        """
        # Reindex HVG row into OS panel (zero-pad unmapped genes)
        x_os = np.zeros(n_os, dtype=np.float32)
        valid = hvg_to_os_arr >= 0
        x_os[hvg_to_os_arr[valid]] = x_row[valid]
        # Pre-normalize='F': log1p(normalize(x))
        total = float(x_os.sum())
        if total <= 0:
            return np.zeros(4 * 768, dtype=np.float32)
        x_norm = np.log1p(x_os / total * 1e4)
        # Append [tgthighres=4.0, log10(total)] (scFoundation t4 default)
        feat = np.concatenate([x_norm, [4.0, np.log10(max(total, 1.0))]])
        pretrain_gene_x = torch.tensor(
            feat, dtype=torch.float32, device="cuda"
        ).unsqueeze(0)
        data_gene_ids = torch.arange(
            19266, device=pretrain_gene_x.device
        ).repeat(pretrain_gene_x.shape[0], 1)
        value_labels = pretrain_gene_x > 0
        x, x_padding = gatherData(pretrain_gene_x, value_labels,
                                   pad_token_id)
        position_gene_ids, _ = gatherData(
            data_gene_ids, value_labels, pad_token_id
        )
        with torch.no_grad():
            x = model_obj.token_emb(
                torch.unsqueeze(x, 2).float(), output_weight=0
            )
            position_emb = model_obj.pos_emb(position_gene_ids)
            x = x + position_emb
            geneemb = model_obj.encoder(x, x_padding)
            geneemb1 = geneemb[:, -1, :]
            geneemb2 = geneemb[:, -2, :]
            geneemb3, _ = torch.max(geneemb[:, :-2, :], dim=1)
            geneemb4 = torch.mean(geneemb[:, :-2, :], dim=1)
            merged = torch.cat([geneemb1, geneemb2, geneemb3, geneemb4],
                               dim=1)
        return merged.squeeze(0).float().cpu().numpy()

    def _embed_batch(X_batch: np.ndarray) -> np.ndarray:
        out = np.stack([_embed_one(X_batch[i]) for i in range(X_batch.shape[0])])
        return out

    BATCH = 4
    probe = _embed_batch(X_all[:1])
    H = int(probe.shape[1])
    embs = np.zeros((n_cells, H), dtype=np.float32)
    embs[0] = probe[0]
    t_emb = time.time()
    for i in range(1, n_cells, BATCH):
        embs[i:i+BATCH] = _embed_batch(X_all[i:i+BATCH])
        if ((i - 1) // BATCH) % 100 == 0:
            print(f"[scfoundation confpert] embedded {i+BATCH}/{n_cells} "
                  f"({time.time()-t_emb:.0f}s)", flush=True)
    print(f"[scfoundation confpert] all embeddings done in "
          f"{time.time()-t_emb:.0f}s (hidden={H})", flush=True)

    # Train/test split over perturbations
    pert_keys = list(ds.X_pert.keys())
    rng = np.random.RandomState(42)
    rng.shuffle(pert_keys)
    n_split = max(2, len(pert_keys) // 2)
    calib_keys = pert_keys[:n_split]
    test_keys = pert_keys[n_split:]
    train_perts = set(calib_keys)

    labels = adata.obs[pert_col].astype(str).values
    train_mask = np.array(
        [(L in train_perts) or is_ctrl_mask[i]
         for i, L in enumerate(labels)], dtype=bool
    )
    train_idx = np.where(train_mask)[0]
    test_idx = np.where(~train_mask)[0]
    print(f"[scfoundation confpert] train cells={train_idx.size} "
          f"test cells={test_idx.size}", flush=True)

    ridge = Ridge(alpha=1.0, fit_intercept=True)
    ridge.fit(embs[train_idx], X_all[train_idx])
    print(f"[scfoundation confpert] ridge fit; coef={ridge.coef_.shape}",
          flush=True)

    train_var = safe_train_variance(
        X_all, train_indices=train_idx, test_indices=test_idx, min_train_n=30
    )

    pred_pops, obs_pops = {}, {}
    t0 = time.time()
    for pi, test_pert in enumerate(test_keys):
        mask = labels == test_pert
        n_obs = int(mask.sum())
        if n_obs < 10:
            continue
        X_obs = X_all[mask]
        emb_test_pert = embs[mask].mean(axis=0)
        mu_pred = ridge.predict(emb_test_pert[None, :])[0].astype(np.float64)
        samples = sample_from_train_distribution(
            mu_per_pert=mu_pred, train_variance=train_var,
            n_cells=n_obs, rng=np.random.default_rng(42 + pi),
        ).astype(np.float32)
        pred_pops[test_pert] = samples
        obs_pops[test_pert] = X_obs.astype(np.float32)
    print(f"[scfoundation confpert] {len(pred_pops)} test perts evaluated",
          flush=True)

    out = _run_conformal_eval(
        pred_pops, obs_pops, dataset, "scfoundation",
        time.time() - t0,
        extra_config={"model_repo": model_repo,
                      "embedding_baseline_mode": True,
                      "os_gene_coverage": int(n_mapped)},
    )
    volume.commit()
    print(f"[scfoundation confpert] done.", flush=True)
    return out


@app.local_entrypoint()
def scfoundation_calibrate(dataset: str = "norman",
                            model_repo: str = "genbio-ai/scFoundation"):
    """scFoundation + ConfPert calibration via .spawn() under --detach."""
    fc = _scfoundation_calibrate_fn.spawn(dataset=dataset, model_repo=model_repo)
    print(f"Spawned scFoundation calibrate: {fc.object_id}")


# ---------------------------------------------------------------------------
# CellFM (Zeng et al. 2024, biomed-AI/CellFM)
# 800M-param RetNet on MindSpore. Largest open scFM. Stress-tests H1 capacity.
# ---------------------------------------------------------------------------

CELLFM_IMAGE = _add_local(
    _BASE_DEPS.apt_install("git")
              .env({
                  "MPLBACKEND": "Agg",
                  "HF_HUB_ENABLE_HF_TRANSFER": "1",
              })
              .pip_install(
                  # CellFM serving is MindSpore-native. The HF mirror provides
                  # only weights; inference still requires the upstream code.
                  # MindSpore install on Modal A100 is the known integration
                  # risk (see predictors_v2_stubs.CELLFM.notes). If install
                  # breaks, fall back to a HF-mirrored PyTorch port should it
                  # exist by Phase 2B.2 window.
                  "transformers>=4.40,<4.46",
                  "huggingface_hub>=0.20,<1.0",
                  "hf_transfer>=0.1.6",
                  # mindspore is NOT pinned at image build time — the import
                  # in _cellfm_calibrate_fn will raise a clear error directing
                  # the user to the MindSpore install matrix for CUDA 12.x.
              )
)


@app.function(
    cpu=4.0,
    memory=24576,
    timeout=60 * 60 * 1,
    volumes={VOLUME_MOUNT: volume},
    image=CELLFM_IMAGE,
    secrets=[modal.Secret.from_name("huggingface")],
)
def _cellfm_prestage_fn(model_repo: str = "biomed-AI/cellfm-800M"):
    """Pre-stage CellFM 800M checkpoint via HF snapshot_download.

    Cost: ~$0.50 CPU (large download).
    """
    print(f"[cellfm-prestage] FUNCTION ENTRY model_repo={model_repo}", flush=True)
    import os
    from huggingface_hub import snapshot_download

    repo_basename = model_repo.split("/")[-1]
    dest_dir = f"{VOLUME_MOUNT}/cellfm_ckpt/{repo_basename}"
    sentinel = f"{dest_dir}/.prestage_complete"

    if os.path.exists(sentinel):
        return {"dest_dir": dest_dir, "skipped": True}

    hf_token = os.environ.get("HF_TOKEN") or os.environ.get("HUGGING_FACE_HUB_TOKEN")
    os.makedirs(dest_dir, exist_ok=True)
    try:
        snapshot_download(repo_id=model_repo, token=hf_token,
                          local_dir=dest_dir, max_workers=1)
    except Exception as exc:
        raise RuntimeError(
            f"snapshot_download({model_repo}) failed: {exc!r}. "
            f"CellFM HF mirror may require license click-through at "
            f"https://huggingface.co/{model_repo}."
        ) from exc

    with open(sentinel, "w") as fh:
        fh.write("prestage complete\n")
    volume.commit()
    return {"dest_dir": dest_dir, "skipped": False}


@app.local_entrypoint()
def cellfm_prestage(model_repo: str = "biomed-AI/cellfm-800M"):
    """Pre-stage CellFM checkpoint."""
    fc = _cellfm_prestage_fn.spawn(model_repo=model_repo)
    print(f"Spawned CellFM pre-stage: {fc.object_id}")


@app.function(
    cpu=4.0,
    memory=32768,
    gpu="A100",
    timeout=60 * 60 * 4,
    volumes={VOLUME_MOUNT: volume},
    image=CELLFM_IMAGE,
    secrets=[modal.Secret.from_name("huggingface")],
)
def _cellfm_calibrate_fn(dataset: str = "norman",
                          n_eval_perts: int = 30,
                          n_top_perturbations: int = 50,
                          model_repo: str = "biomed-AI/cellfm-800M"):
    """CellFM perturbation prediction + ConfPert calibration (STUB).

    Phase 2B.2 implementation flow:
      1. Confirm MindSpore (or PyTorch port) installs successfully in the
         CELLFM_IMAGE. If MindSpore breaks on Modal A100, defer per
         predictors_v2_stubs.CELLFM.notes.
      2. Load 800M checkpoint, extract per-cell embeddings.
      3. Train perturbation head on TRAIN labels.
      4. Sample via safe_train_variance + sample_from_train_distribution.
      5. Run `_run_conformal_eval`.

    Heavyweight predictor cap: $80/(predictor, dataset) per pre-reg §14.5
    (highest of the 4 heavyweight predictors due to 800M param count).
    """
    raise NotImplementedError(
        "CellFM wrapper is a Phase 2B.2 deliverable. See inline TODO. "
        "Per pre-reg v2 §1.6: report as N/A until wrapper passes user-supervised "
        "GPU smoke test. Use confpert.heavyweight_helpers.safe_train_variance."
    )


@app.local_entrypoint()
def cellfm_calibrate(dataset: str = "norman",
                      model_repo: str = "biomed-AI/cellfm-800M"):
    """CellFM + ConfPert calibration via .spawn() under --detach."""
    fc = _cellfm_calibrate_fn.spawn(dataset=dataset, model_repo=model_repo)
    print(f"Spawned CellFM calibrate: {fc.object_id}")
