# Circularity Detection Score (CDS)

Detect circular dependencies between benchmark endpoints and prediction features in transcriptomic studies.

## Quick Start

```bash
pip install .                                              # from a checkout of this package
# or: pip install git+https://github.com/huangzhixian324-dev/cds.git@v1.2.1
cds --demo
```

## What is CDS?

CDS quantifies how much a gene-expression-derived response endpoint overlaps with the features available for prediction. This is critical for benchmark studies: if your "response" label is itself defined by gene expression (e.g., cytolytic activity = GZMA + PRF1), then ML methods will trivially learn to predict it, achieving inflated AUROC values that do not transfer to clinical endpoints.

CDS ranges from 0 (no circularity) to 1 (maximum circularity):
- **CDS > 0.70**: HIGH risk — use a clinical endpoint instead
- **CDS 0.30-0.70**: MODERATE risk — interpret results cautiously
- **CDS < 0.30**: LOW risk — endpoint likely independent

## Usage

### Command Line

```bash
# Demo with synthetic data
cds --demo

# Real data
cds --expr expression.tsv --labels response.csv --endpoint-genes GZMA,PRF1

# Generate HTML report
cds --expr data.tsv --labels labels.csv --endpoint-genes GZMA,PRF1 --html report.html
```

### Python API

```python
from cds import circularity_detection_score
import numpy as np

# Simulated data with GZMA/PRF1-defined endpoint
X = np.random.lognormal(4, 1, (100, 5000))
y = (X[:,0] + X[:,1] > np.median(X[:,0] + X[:,1])).astype(int)
genes = ["GZMA","PRF1"] + [f"GENE_{i}" for i in range(4998)]

result = circularity_detection_score(X, y, ["GZMA","PRF1"], genes)
print(f"CDS = {result['CDS']:.3f} [{result['risk_level']}]")
# Output: CDS = 0.852 [HIGH]
```

## Reference

If you use CDS in your research, please cite:

> SPATBench: Methodological Choices Can Outweigh Algorithmic Choices in Transcriptomic Immunotherapy Prediction Benchmarking at Current Sample Sizes. (2026)

## License

MIT
