# PFCI: Penalized Fast Causal Inference for High-Dimensional Structure Learning

PFCI implements **Penalized Fast Causal Inference (PFCI)**, a scalable two-stage procedure for learning graphical structures in high-dimensional settings with potential latent variables and selection bias.

The method combines:

- **Graphical lasso screening** to obtain a sparse super-skeleton  
- **Constrained Fast Causal Inference (FCI)** for orientation and refinement  

This enables computationally efficient structure learning while preserving theoretical guarantees under sparsity assumptions.

---

## 📦 Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
install.packages("BiocManager")
BiocManager::install(c("graph", "RBGL", "ggm", "pcalg"))
devtools::install_github("SamhitaPal3/PFCI")

library(PFCI)

sim <- simulate_pfci_toy()
fit <- pfci_fit(sim$X)
met <- pfci_metrics(sim, fit)
met

BiocManager::install("Rgraphviz")
plot_pag(fit)
```
# Reference
Pal, S., Ghosh, D., & Yang, S. (2025). Penalized FCI for Causal Structure Learning in a Sparse DAG for Biomarker Discovery in Parkinson's Disease. arXiv preprint arXiv:2507.00173.
