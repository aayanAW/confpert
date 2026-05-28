FROM python:3.11-slim

WORKDIR /submission

# Bake in the minimum runtime; participants extend with their own predictor
# weights / framework deps (torch, jax, scgpt, etc.)
RUN pip install --no-cache-dir \
    numpy>=1.26 \
    scipy>=1.11 \
    scanpy>=1.10

COPY predictor_template.py /submission/predictor.py

# Submissions OVERRIDE predictor.py with their own implementation of
# predict_samples(X_ctrl_test, perturbation, n_cells) -> np.ndarray.
# The Codabench scorer mounts /input + /output and runs:
#   python -m predictor

ENTRYPOINT ["python", "-m", "predictor"]
