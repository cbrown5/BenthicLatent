#' Plot flow and distance effects
#' for multinomial poisson model
#'
#' @details
#' This has been cut from the paper now, b/c I decided not to present the
#' cumulative distance results.  
#' @author Christopher J. Brown
#' @rdname plot_effects
#' @export

plot_effects <- function(fmindist, fmindist2, fcumdist){
	load(fmindist)
	mcoutmd <- mcout3
	load(fmindist2)
	mcoutmd2 <- mcout3
	load(fcumdist)


	gamma_md <- jagstools::jagsresults(x = mcoutmd, params = c('gamma'))
	gammaf_md <- jagstools::jagsresults(x = mcoutmd, params = c('gammaf'))

	gamma_md2 <- jagstools::jagsresults(x = mcoutmd2, params = c('gamma'))
	gammaf_md2 <- jagstools::jagsresults(x = mcoutmd2, params = c('gammaf'))

	gamma_cd <- jagstools::jagsresults(x = mcout3, params = c('gamma'))
	gammaf_cd <- jagstools::jagsresults(x = mcout3, params = c('gammaf'))

	#yvals <- c(1,2,3,4,6,7)
	yvals <- c(1,2,3,4)
	#labs <- c('Min. \n distance 6 LVs', 'High flow \n (min. distance) 6 LVs',
	#'Min. \n distance 2 LVs', 'High flow \n (min. distance) 2 LVs',
#'Cumulative \n distance','High flow  \n (cumulative \n distance)')
	labs <- c('Min. \n distance 2 LVs', 'High flow \n (min. distance) 2 LVs',
	'Cumulative \n distance','High flow  \n (cumulative \n distance)')

	#xvals <- c(gamma_md2[1], gammaf_md2[1],gamma_md[1], gammaf_md[1], gamma_cd[1], gammaf_cd[1])
	#xlwr <- c(gamma_md2[3], gammaf_md2[3],gamma_md[3], gammaf_md[3], gamma_cd[3], gammaf_cd[3])
	#xupr <- c(gamma_md2[7], gammaf_md2[7],gamma_md[7], gammaf_md[7], gamma_cd[7], gammaf_cd[7])
	xvals <- c(gamma_md[1], gammaf_md[1], gamma_cd[1], gammaf_cd[1])
	xlwr <- c(gamma_md[3], gammaf_md[3], gamma_cd[3], gammaf_cd[3])
	xupr <- c(gamma_md[7], gammaf_md[7], gamma_cd[7], gammaf_cd[7])

	par(mar = c(5,7,3,2))
	plot(xvals, yvals, xlim = c(min(xlwr)*1.5, max(xupr)*1.5),
 yaxt = 'n', bty = 'n',
  ylab = '', xlab = 'Effect size', pch = 16, col = grey(0.5, 1),
	 cex =0.7, xaxs = 'i')
	arrows(xlwr, yvals, xupr, yvals, len = 0,
 col = grey(0.5, 0.5), lwd = 2)
 abline(v=0, lty = 2)
	axis(2, at = yvals, labels = labs, las=1, cex.axis = 0.6)
}
