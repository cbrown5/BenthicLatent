#' Plot loadings on latent variables 
#' for multinomial model
#' 
#' @author Christopher J. Brown
#' @rdname plot_latent_mnm
#' @export

plot_latent_mnm <- function(smc,  habnams, nlv,wqlv = 1, xlim = c(-3, 3), xlab = '', ylab = ''){
	dimnam <- attr(smc$stat, 'dimnames')[[1]]
	ihab <- grep('all.params',dimnam)
	vhab <- smc$statistics[ihab,1]
	names(vhab) <- paste0(rep(habnams, nlv), 
	rep(c(paste0(rep('_lv',nlv), 1:nlv)), each =length(ihab)/(nlv)))

	iwq <- grep(paste0('lv',wqlv),names(vhab))

	sppnums <- 1:(length(ihab)/(nlv))
	
	iord <- order(vhab[iwq])
	ilord <-iwq[iord]
	
	plot(vhab[ilord], sppnums, yaxt = 'n', bty = 'n',
	 xlim = xlim, ylab = ylab, xlab = xlab, pch = 16, col = grey(0.5, 1), 
	 cex =0.7)
	abline(v=0)
	arrows(x0 = smc$quantiles[ihab[ilord],1], y0=sppnums,
	 x1 = smc$quantiles[ihab[ilord],5], y1=sppnums, len = 0, col = grey(0.5, 0.5), lwd = 2)
	axis(2, at = sppnums, labels = habnams[iord], las=1, cex.axis = 0.5)
}
