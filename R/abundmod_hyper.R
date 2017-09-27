#' Predict abundances from hyper-parameters
#'
#' @author Christopher J. Brown
#' @rdname abundmod_hyper
#' @export


abundmod_hyper <- function(i, nlvs, n, nhabs, X1, X2, alpha, beta, int, gamma1, gamma2){
	amult <- rep(0,nhabs)
	for (j in 1:nlvs){
		nu <- 0#rnorm(1,0,1)
		amult <- amult +
		(alpha[((j * nhabs) - nhabs + 1):(j * nhabs)] * nu)
		}
	nuwq <- (gamma1 * X1[i]) + (gamma2 * X2[i])# + rnorm(1,0,1)

	exp(beta*nuwq + amult + int)
	}
