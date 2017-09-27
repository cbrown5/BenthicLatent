#' Predict abundances from latent variable estimates
#'
#' @author Christopher J. Brown
#' @rdname abundmod
#' @export


abundmod <- function(i, nlvs, n, nhabs, nuwq, nu, alpha, beta, int){
	amult <- rep(0,nhabs)
	for (j in 1:nlvs){
		amult <- amult +
		(alpha[((j * nhabs) - nhabs + 1):(j * nhabs)] *
		nu[(j*n - n +i)])
		}
	exp(beta*nuwq[i] + amult + int)
	}
