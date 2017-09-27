#' Plot loadings on fixed latent variable
#' for multinomial poisson model
#'
#' @author Christopher J. Brown
#' @rdname plot_beta_mnm
#' @export

plot_beta_mnm <- function(smc,  habnams, xlim = NULL, xlab = '', ylab = ''){

	dimnam <- attr(smc$stat, 'dimnames')[[1]]
	ihab <- grep('beta',dimnam)
	vhab <- smc$statistics[ihab,1]
	names(vhab) <- habnams

	sppnums <- 1:length(habnams)
	iord <- order(vhab)

	if (is.null(xlim)){
		xmin <- max(abs(range(c(smc$quantiles[ihab,1], smc$quantiles[ihab,5]))))
		xlim <- c(-xmin, xmin)
		}

	plot(vhab[iord], sppnums, yaxt = 'n', bty = 'n',
	 xlim = xlim, ylab = ylab, xlab = xlab, pch = 16, col = grey(0.5, 1),
	 cex =0.7)
	abline(v=0)
	arrows(x0 = smc$quantiles[ihab[iord],1], y0=sppnums,
	 x1 = smc$quantiles[ihab[iord],5], y1=sppnums, len = 0, col = grey(0.5, 0.5), lwd = 2)
	xy <- axis(2, at = sppnums, labels = NA, las=1)
	text(xlim[1]*1.1, xy, habnams[iord], srt = 35, cex = 0.8, xpd = NA, pos=2)
	return(iord)
}
