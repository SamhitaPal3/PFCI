# PFCI

Penalized Fast Causal Inference for high-dimensional structure learning.

## Installation
devtools::install_github("SamhitaPal3/PFCI")

## Basic Workflow
sim <- simulate_pfci_toy(...)
fit <- pfci_fit(sim$X)
metrics(sim, fit)

## Simulation and Reproduction Scripts
See inst/scripts/
