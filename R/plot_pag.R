#' Plot a PAG returned by PFCI
#'
#' @param fit A pfci_fit object returned by pfci_fit(), or a pcalg fci object.
#' @param ... Additional arguments passed to pcalg plotting.
#'
#' @export
plot_pag <- function(fit, ...) {

  if (inherits(fit, "pfci_fit")) {
    pag_obj <- fit$pag
  } else {
    pag_obj <- fit
  }

  # Explicitly call pcalg plot method
  pcalg::plot(pag_obj, ...)
}
