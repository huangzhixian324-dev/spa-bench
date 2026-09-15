FROM python:3.11-slim

LABEL org.opencontainers.image.title="SPATBench"
LABEL org.opencontainers.image.description="Benchmark framework isolating methodological design choices (response definition, feature selection, endpoint circularity) in transcriptomic immunotherapy prediction"
LABEL org.opencontainers.image.version="v33"
LABEL org.opencontainers.image.licenses="MIT"

# System dependencies: h5py needs the HDF5 headers at build time
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Python dependencies (core analysis stack; see requirements.txt)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Project code. Data and outputs are excluded via .dockerignore; the archived
# inputs live in the Zenodo record and can be mounted at runtime.
COPY . .

# Default action: the end-to-end synthetic self-test (exit code 0 expected).
# Reproducing the reported tables additionally requires the archived inputs:
#   docker run -v $PWD/results:/app/results spa-bench:v33 \
#       python -W ignore scripts/rerun_v33.py
CMD ["python", "-W", "ignore", "scripts/validate.py"]
