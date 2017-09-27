#' Fit new LV model
#' 
#' @author Christopher J. Brown
#' @rdname fit_lvmod
#' @export


fit_lvmod <- function(lv_input, dmetric, num.lv = 2, 
nburnin = 5000, niter = 100000, nthin = 10){

	load.module("glm")
	setwd('~/Code/BenthicLatent')
		
	#mindist model
	jags <- jags.model('bug/lv_multinom_dists.bug', 
	data = list('n' = lv_input$N, 
	'p' =  lv_input$p,
	'num.lv' = 2, 
	'd' =  dmetric, 
	'y' =  lv_input$y, 
	'flow' =  lv_input$flow),
	n.chains = 1)
		
	#Burn in
	update(jags, nburnin)
	
	#Extract samples
	niter <- niter
	nthin <- nthin
	nsamp <- round(niter/nthin)
	lvmod_out <- coda.samples(jags, variable.names=c("alpha","alphaf","lvs","int", "all.params"), n.iter=niter, thin = nthin)
	return(lvmod_out)
	}