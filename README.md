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
devtools::install_github("SamhitaPal3/PFCI")

library(PFCI)

sim <- simulate_pfci_toy()
fit <- pfci_fit(sim$X)
met <- pfci_metrics(sim, fit)
met
```
