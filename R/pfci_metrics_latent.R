#' Metrics for latent simulation using oracle-FCI truth (skeleton only)
#'
#' Designed for the 3-line workflow:
#'   sim <- simulate_with_latent(...)
#'   fit <- pfci_fit(sim$X, ...)
#'   met <- metrics_with_latent(sim, fit)
#'
#' Returns only: SHD, F1_total, MCC, Time.
#'
#' @param sim Output from simulate_with_latent().
#' @param fit Output from pfci_fit() (must contain $amat and $time$total).
#'
#' @return A named list with SHD, F1_total, MCC, Time.
#' @export
metrics_with_latent <- function(sim, fit) {
  stopifnot(is.list(sim), is.list(sim$truth), is.matrix(sim$truth$skel))
  stopifnot(is.list(fit), !is.null(fit$amat))

  true_skel <- sim$truth$skel
  est_skel  <- .amat_to_skeleton_latent(as(fit$amat, "matrix"))

  SHD <- .skel_shd_latent(true_skel, est_skel)
  fm  <- .skel_f1_mcc_latent(true_skel, est_skel)

  list(
    SHD = as.numeric(SHD),
    F1_total = unname(fm["F1"]),
    MCC = unname(fm["MCC"]),
    Time = if (!is.null(fit$time$total)) as.numeric(fit$time$total) else NA_real_
  )
}

# -----------------------
# INTERNAL HELPERS (NOT EXPORTED)
# -----------------------

.amat_to_skeleton_latent <- function(amat) {
  A <- (amat != 0L) * 1L
  A <- ((A + t(A)) > 0L) * 1L
  diag(A) <- 0L
  A
}

.skel_shd_latent <- function(true_skel, est_skel) {
  p <- nrow(true_skel)
  d <- 0L
  for (i in 1L:(p - 1L)) for (j in (i + 1L):p) {
    if (true_skel[i, j] != est_skel[i, j]) d <- d + 1L
  }
  d
}

.skel_f1_mcc_latent <- function(true_skel, est_skel) {
  ut <- upper.tri(true_skel, diag = FALSE)
  tU <- true_skel[ut]
  eU <- est_skel[ut]

  # Force numeric doubles (prevents integer overflow)
  TP <- sum(tU == 1 & eU == 1, na.rm = TRUE) * 1.0
  FP <- sum(tU == 0 & eU == 1, na.rm = TRUE) * 1.0
  FN <- sum(tU == 1 & eU == 0, na.rm = TRUE) * 1.0
  TN <- sum(tU == 0 & eU == 0, na.rm = TRUE) * 1.0

  prec <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  rec  <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  f1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0

  denom <- (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
  mcc <- if (is.finite(denom) && denom > 0) (TP * TN - FP * FN) / sqrt(denom) else 0

  c(F1 = f1, MCC = mcc)
}
