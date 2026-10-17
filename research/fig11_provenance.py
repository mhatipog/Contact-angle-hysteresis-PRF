"""Figure 11 archived bisection lattice audit; no physical solver invoked."""
from pathlib import Path
import numpy as np
from scipy.io import loadmat

def audit(path=Path("Data/numerics_comparison_data.mat")):
    source = loadmat(path)
    V = float(source["V"].item())
    table = np.asarray(source["xu0_escape"], dtype=float)
    finite = table[np.isfinite(table)]
    if table.shape != (5, 11) or finite.size != 41:
        raise AssertionError("Unexpected archived Figure 11 dataset")
    output = {}
    for name, lo, hi in [("committed_generator", V+.01, 1.0),
                         ("paper_text", V+.03, .97)]:
        iterations = int(np.ceil(np.log2((hi-lo)/.005)))
        step = (hi-lo)/(2**iterations)
        distances = np.abs(((finite-lo)/step)-np.round((finite-lo)/step))*step
        output[name] = {"iterations": iterations, "step": step,
                        "exact_count": int(np.sum(distances < 1e-10)),
                        "median_error": float(np.median(distances)),
                        "max_error": float(np.max(distances))}
    return output

if __name__ == "__main__":
    results = audit()
    for label, metrics in results.items():
        print(label, metrics)
    assert results["committed_generator"]["exact_count"] == 41
    assert results["paper_text"]["exact_count"] == 0
