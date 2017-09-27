#' Predict abundances from hyper parameters across all sites
#'
#' @author Christopher J. Brown
#' @rdname mapabund_hyper
#' @export


mapabund_hyper <- function(nlvs, n, nhabs, X1, X2, alpha, beta, int, gamma1, gamma2){
	purrr::map(1:n, ~abundmod_hyper(.x,nlvs, n, nhabs, X1, X2, alpha, beta, int, gamma1, gamma2))
}
