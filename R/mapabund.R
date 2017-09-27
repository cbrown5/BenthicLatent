#' Predict abundances from latent variable estimates across all sites
#'
#' @author Christopher J. Brown
#' @rdname abundmod
#' @export


mapabund <- function(nlvs, n, nhabs, nuwq, nu, alpha, beta, int){
	map(1:n, ~abundmod(.x,nlvs, n, nhabs, nuwq, nu, alpha, beta, int))
}
