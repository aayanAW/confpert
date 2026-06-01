"""Modal CPU runner for the weighted-CP cross-cell-line experiment (Track B Task 5 Step 8).

Runs scripts/run_weighted_cp_xcl.py's logic on a Modal CPU worker with the
causeflow-artifacts volume mounted, so the large Replogle h5ads never touch the
local disk (the local run filled it). Reuses the committed, smoke-tested `run()`
verbatim -- this wrapper only mounts the data and ships the result back.

Launch:
  modal run --detach scripts/modal_weighted_cp_xcl.py::weighted_cp_xcl
Then locally write the returned dict to results/weighted_cp_xcl.json and commit.
"""
from __future__ import annotations

from pathlib import Path

import modal

APP_NAME = "confpert-wcp"

IMAGE = (
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
    .add_local_dir(str(Path(__file__).resolve().parent.parent / "src"),
                   remote_path="/root/src")
    .add_local_dir(str(Path(__file__).resolve().parent),
                   remote_path="/root/scripts")
)

app = modal.App(APP_NAME, image=IMAGE)
volume = modal.Volume.from_name("causeflow-artifacts", create_if_missing=True)


@app.function(volumes={"/artifacts": volume}, timeout=60 * 60, retries=0)
def weighted_cp_xcl(predictor: str = "ahlmann_bilinear_ridge",
                     noise: str = "C_per_gene_marginal", seed: int = 42) -> dict:
    import os
    import sys

    # The committed run() reads ROOT/data/replogle_*.h5ad where ROOT=/root (the
    # script lives at /root/scripts). Point /root/data at the mounted volume's
    # /artifacts/data so no download is needed.
    os.makedirs("/root", exist_ok=True)
    if not os.path.exists("/root/data"):
        os.symlink("/artifacts/data", "/root/data")
    print("[wcp-modal] /root/data ->", os.path.realpath("/root/data"), flush=True)
    print("[wcp-modal] data files:",
          [f for f in os.listdir("/root/data") if "replogle" in f], flush=True)

    sys.path.insert(0, "/root/src")
    sys.path.insert(0, "/root/scripts")
    from run_weighted_cp_xcl import run  # committed, smoke-tested logic

    res = run(predictor, noise, seed=seed)

    # Persist a copy on the volume too (belt-and-suspenders; the return value is
    # the authoritative channel back to the launcher).
    import json
    os.makedirs("/artifacts/results", exist_ok=True)
    with open("/artifacts/results/weighted_cp_xcl.json", "w") as f:
        json.dump(res, f, indent=2)
    volume.commit()
    return res


@app.local_entrypoint()
def main(predictor: str = "ahlmann_bilinear_ridge",
         noise: str = "C_per_gene_marginal", seed: int = 42):
    import json

    res = weighted_cp_xcl.remote(predictor=predictor, noise=noise, seed=seed)
    out = Path(__file__).resolve().parent.parent / "results" / "weighted_cp_xcl.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2))
    print(f"[wcp-modal] wrote {out}")
    print(f"[wcp-modal] ess={res.get('weight_ess'):.2f} "
          f"n_perts={res.get('n_common_perts')}")
    for sn, rec in res.get("per_score", {}).items():
        a = rec["alpha_0.2"]
        print(f"  {sn:22s} unw={a['unweighted_coverage']:.3f} "
              f"wtd={a['weighted_coverage']:.3f}")
