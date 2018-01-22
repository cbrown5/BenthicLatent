#' Plot loadings on latent variables
#' for multinomial poisson model
#'
#' @author Christopher J. Brown
#' @rdname plot_alpha_mnm
#' @export

plot_alpha_mnm <- function(smc,  habnams, lv=1,nlv, iord, xlim = NULL, xlab = '', ylab = ''){
	dimnam <- attr(smc$stat, 'dimnames')[[1]]
	ihab <- grep('alpha',dimnam)
	vhab <- smc$statistics[ihab,1]
	names(vhab) <- paste0(rep(habnams, nlv),
	rep(c(paste0(rep('_lv',nlv), 1:nlv)), each =length(ihab)/(nlv)))

	sppnums <- 1:(length(ihab)/(nlv))

	iplot <- grep(lv, names(vhab))
	#iord <- order(vhab[iplot])

	if (is.null(xlim)){
		xmin <- max(abs(range(c(smc$quantiles[ihab[iplot][iord],1], smc$quantiles[ihab[iplot][iord],5]))))
		xlim <- c(-xmin, xmin)
		}

	plot(vhab[iplot][iord], sppnums, yaxt = 'n', bty = 'n',
	 xlim = xlim, ylab = ylab, xlab = xlab, pch = 16, col = "black",
	 cex =0.9)
	abline(v=0)
	arrows(x0 = smc$quantiles[ihab[iplot][iord],1], y0=sppnums,
	 x1 = smc$quantiles[ihab[iplot][iord],5], y1=sppnums, len = 0, col = grey(0, 0.5), lwd = 1.5)
	#axis(2, at = sppnums, labels =NA, las=1, cex.axis = 0.5)
	#axis(2, at = sppnums, labels = habnams[iord], las=1, cex.axis = 0.5)
}
