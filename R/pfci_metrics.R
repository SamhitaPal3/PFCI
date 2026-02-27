#' Compute PFCI metrics from a simulation object and a pfci_fit output
#'
#' Designed for the 3-line workflow:
#'   sim <- simulate_pfci_toy(...)
#'   fit <- pfci_fun(sim$X, ...)
#'   met <- pfci_metrics(sim, fit)
#'
#' Default metrics compare estimated PAG adjacency (skeleton) to the generating DAG skeleton.
#'
#' If compute_marks=TRUE and sim$truth$amat exists, it also reports mark-level F1s:
#'  - F1_dir  (->)
#'  - F1_oDir (o->)
#'  - F1_bidir (<->)
#'  - F1_circ (o-o)
#'  - F1_arrow (arrowheads)
#'  - F1_tail  (tails)
#'
#' @param sim Output from simulate_pfci_toy().
#' @param fit Output from pfci_fun()/pfci_fit() with at least $amat and $time$total.
#' @param compute_marks Logical. If TRUE, also computes mark-level F1 when truth amat is present.
#'
#' @return A named list of metrics.
#' @export
pfci_metrics <- function(sim, fit, compute_marks = FALSE) {
  stopifnot(is.list(sim), is.matrix(sim$X))
  stopifnot(is.list(fit), !is.null(fit$amat))

  true_skel <- sim$truth$skel
  if (is.null(true_skel) || !is.matrix(true_skel)) {
    stop("sim$truth$skel is missing. Use simulate_pfci_toy() output.")
  }

  est_skel <- .amat_to_skeleton(as(fit$amat, "matrix"))

  shd <- .skel_shd(true_skel, est_skel)
  fm  <- .skel_f1_mcc(true_skel, est_skel)

  out <- list(
    SHD = as.numeric(shd),
    F1_total = unname(fm["F1"]),
    MCC = unname(fm["MCC"]),
    Precision = unname(fm["Precision"]),
    Recall = unname(fm["Recall"]),
    TP = unname(fm["TP"]), FP = unname(fm["FP"]), FN = unname(fm["FN"]), TN = unname(fm["TN"]),
    Time = if (!is.null(fit$time$total)) as.numeric(fit$time$total) else NA_real_,
    time_glasso = if (!is.null(fit$time$glasso)) as.numeric(fit$time$glasso) else NA_real_,
    time_fci = if (!is.null(fit$time$fci)) as.numeric(fit$time$fci) else NA_real_,
    rho = if (!is.null(fit$rho)) as.numeric(fit$rho) else NA_real_
  )

  if (isTRUE(compute_marks)) {
    truth_amat <- sim$truth$amat
    if (!is.null(truth_amat) && is.matrix(truth_amat)) {
      est_amat <- as(fit$amat, "matrix")
      out$F1_dir   <- .class_f1(truth_amat, est_amat, "dir")
      out$F1_oDir  <- .class_f1(truth_amat, est_amat, "oDir")
      out$F1_bidir <- .class_f1(truth_amat, est_amat, "bidir")
      out$F1_circ  <- .class_f1(truth_amat, est_amat, "circ")
      out$F1_arrow <- .endpoint_f1(truth_amat, est_amat, "arrow")
      out$F1_tail  <- .endpoint_f1(truth_amat, est_amat, "tail")
    } else {
      out$F1_dir <- out$F1_oDir <- out$F1_bidir <- out$F1_circ <- NA_real_
      out$F1_arrow <- out$F1_tail <- NA_real_
    }
  }

  out
}

# -------------------------------------------------------------------------
# INTERNAL helpers (NOT exported)
# -------------------------------------------------------------------------

.amat_to_skeleton <- function(amat) {
  A <- (amat != 0L) * 1L
  A <- ((A + t(A)) > 0L) * 1L
  diag(A) <- 0L
  A
}

.skel_shd <- function(true_skel, est_skel) {
  p <- nrow(true_skel)
  d <- 0L
  for (i in 1:(p - 1L)) for (j in (i + 1L):p) {
    if (true_skel[i, j] != est_skel[i, j]) d <- d + 1L
  }
  d
}

# .skel_f1_mcc <- function(true_skel, est_skel) {
#   # count each undirected pair once (upper triangle)
#   ut <- upper.tri(true_skel, diag = FALSE)
#   tvec <- true_skel[ut]
#   evec <- est_skel[ut]
#
#   # force numeric (double) to avoid integer overflow
#   TP <- as.numeric(sum(tvec == 1L & evec == 1L))
#   FP <- as.numeric(sum(tvec == 0L & evec == 1L))
#   FN <- as.numeric(sum(tvec == 1L & evec == 0L))
#   TN <- as.numeric(sum(tvec == 0L & evec == 0L))
#
#   prec <- if ((TP + FP) > 0) TP / (TP + FP) else 0
#   rec  <- if ((TP + FN) > 0) TP / (TP + FN) else 0
#   f1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
#
#   denom <- (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
#
#   # guard against weird numeric issues
#   mcc <- if (is.finite(denom) && denom > 0) {
#     (TP * TN - FP * FN) / sqrt(denom)
#   } else {
#     0
#   }
#
#   c(F1 = f1, MCC = mcc, Precision = prec, Recall = rec,
#     TP = TP, FP = FP, FN = FN, TN = TN)
# }

# INTERNAL (not exported)
.skel_f1_mcc <- function(true_skel, est_skel) {
  # Inputs: symmetric 0/1 matrices, diag=0
  stopifnot(is.matrix(true_skel), is.matrix(est_skel))
  stopifnot(all(dim(true_skel) == dim(est_skel)))

  p <- nrow(true_skel)

  # Use upper triangle to avoid double counting and to avoid huge integer ops
  ut <- upper.tri(true_skel, diag = FALSE)

  tU <- true_skel[ut]
  eU <- est_skel[ut]

  # Force numeric early (prevents integer overflow)
  TP <- sum(tU == 1 & eU == 1, na.rm = TRUE) * 1.0
  FP <- sum(tU == 0 & eU == 1, na.rm = TRUE) * 1.0
  FN <- sum(tU == 1 & eU == 0, na.rm = TRUE) * 1.0
  TN <- sum(tU == 0 & eU == 0, na.rm = TRUE) * 1.0

  prec <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  rec  <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  f1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0

  denom <- (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)

  mcc <- 0
  if (is.finite(denom) && denom > 0) {
    mcc <- (TP * TN - FP * FN) / sqrt(denom)
  } else {
    mcc <- 0
  }

  c(
    F1 = f1, MCC = mcc,
    Precision = prec, Recall = rec,
    TP = TP, FP = FP, FN = FN, TN = TN
  )
}

.safe_prf <- function(TP, FP, FN) {
  prec <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  rec  <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  f1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
  c(precision = prec, recall = rec, f1 = f1)
}

.pair_class <- function(amat, i, j) {
  a <- amat[i, j]
  b <- amat[j, i]
  if (a == 0L && b == 0L) return("none")
  if ((a == 2L && b == 3L) || (a == 3L && b == 2L)) return("dir")
  if ((a == 2L && b == 1L) || (a == 1L && b == 2L)) return("oDir")
  if (a == 2L && b == 2L) return("bidir")
  if (a == 1L && b == 1L) return("circ")
  "other"
}

.class_f1 <- function(true_amat, est_amat, cls = c("dir", "oDir", "bidir", "circ")) {
  cls <- match.arg(cls)
  p <- nrow(true_amat)
  TP <- FP <- FN <- 0L
  for (i in 1:(p - 1L)) for (j in (i + 1L):p) {
    t <- (.pair_class(true_amat, i, j) == cls)
    e <- (.pair_class(est_amat , i, j) == cls)
    if (e && t)  TP <- TP + 1L
    if (e && !t) FP <- FP + 1L
    if (!e && t) FN <- FN + 1L
  }
  .safe_prf(TP, FP, FN)[["f1"]]
}

.endpoint_f1 <- function(true_amat, est_amat, mark = c("arrow", "tail")) {
  mark <- match.arg(mark)
  code <- if (mark == "arrow") 2L else 3L
  p <- nrow(true_amat)
  TP <- FP <- FN <- 0L
  for (i in 1:p) for (j in 1:p) if (i != j) {
    t <- (true_amat[i, j] == code)
    e <- (est_amat [i, j] == code)
    if (e && t)  TP <- TP + 1L
    if (e && !t) FP <- FP + 1L
    if (!e && t) FN <- FN + 1L
  }
  .safe_prf(TP, FP, FN)[["f1"]]
}
