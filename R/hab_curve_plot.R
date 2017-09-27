#' Curves with CIs for signficant habitats
#'
#' @author Christopher J. Brown
#' @rdname hab_curve_plot
#' @export


hab_curve_plot <- function(savnam){
	data(lv_input)
	data(bentbl)
	load(savnam)

	#Covariates
	x <- lv_input$mindist
	flowind <- lv_input$flow
	n <- nrow(lv_input$y)
	nhabs <- ncol(lv_input$y)
	habnams <- dimnames(lv_input$y)[[2]]

	colnams <- dimnames(mcout3[[1]])[[2]]

	ialpha <- stringr::str_detect(colnams, 'alpha')
	ibeta <- stringr::str_detect(colnams, 'beta')
	iint <- stringr::str_detect(colnams, 'int')
	igam <- stringr::str_detect(colnams, '\\bgamma\\b')
	igamf <- stringr::str_detect(colnams, 'gammaf')
	inu <- stringr::str_detect(colnams, '\\bnu\\b')

	num_levels <- sum(ialpha)/nhabs

	mcmat <- as.matrix(mcout3[[1]])

	params <- list(alpha = mcmat[,ialpha],
	beta = mcmat[,ibeta],
	int = mcmat[,iint],
	gamma = mcmat[,igam],
	gammaf = mcmat[,igamf],
	nu = mcmat[,inu])

	nsamp <- nrow(mcmat)

	nx <- 20
	xvals <- seq(min(x), max(x), length.out = nx)
	predflow <- rep(0, nx)

	df <- purrr::map(1:nsamp, ~mapabund_hyper(num_levels, nx, nhabs, xvals, predflow, params$alpha[.x,],
		params$beta[.x,], params$int[.x,], params$gamma[.x], params$gammaf[.x])) %>%
		purrr::transpose()
	#so this creates a list structured in order of - predval - sample - habitat.

	get_pred_ci <- function(i, df, nsamp, habnams, quants = c(0.025, 0.5, 0.975)){
		df2 <- df[[i]] %>% unlist() %>% matrix(nrow = nsamp, byrow = T) %>%
		data.frame()
		proppred <- df2/rowSums(df2)
		df3 <- apply(proppred, 2, quantile, probs = quants) %>%
		t() %>% data.frame() %>% cbind(ix = i, hab = habnams)
		df3
	}


	smc <- summary(mcout3)
	betaquant <- smc$quantiles[ibeta,]

	issig <- purrr::map_lgl(1:nhabs, ~sign(betaquant[.x,1]) == sign(betaquant[.x,5]) )
	sighab <- habnams[issig]
	betainc <- betaquant[issig,3]
	datord <- data.frame(hab = habnams[issig], beta = betainc)

	dfci <- purrr::map(1:nx, ~get_pred_ci(.x, df, nsamp, habnams)) %>% dplyr::bind_rows()

	dfci$xvals <- xvals[as.numeric(dfci$ix)]
	dfci$dist <- (dfci$xvals * sd(lv_input$mindists)) + mean(lv_input$mindists)
	dffilt <- dplyr::filter(dfci, hab %in% sighab) %>%
	dplyr::left_join(bentbl, by = c('hab' = 'CODE')) %>%
	dplyr::left_join(datord)

	dffilt$CATEGORY[is.na(dffilt$CATEGORY)] <- dffilt$hab[is.na(dffilt$CATEGORY)]
	dffilt$CATEGORY[dffilt$hab == 'branching_dead'] <- 'Dead branching coral'
	dffilt$CATEGORY[dffilt$hab == 'coral_rare'] <- 'Rare corals'
	dffilt$CATEGORY[dffilt$hab == 'other'] <- 'Other habitats'
	dffilt$CATEGORY[dffilt$hab == 'unconsolidated'] <- 'Soft sediment'
	dffilt$CATEGORY[dffilt$hab == 'coral_submassive'] <- 'Submassive coral'
	dffilt$cat2 <- reorder(dffilt$CATEGORY, dffilt$beta, min)
	#Convert distances back to km
	dffilt$distkm <- ((dffilt$dist * lv_input$sddist) + lv_input$mndist)/1000

	# par(mfrow = c(2,5), mar = c(2,2,2,2), las = 1)
	# for (i in 1:length(sighab)){
		# df <- filter(dffilt, hab == sighab[i])
		# plot(df$dist, df$X50., type = 'l', ylim = c(0, 0.3), xlab = '', ylab = '')
		# addpoly(df$dist, y0 = df$X2.5., y = df$X97.5., col = grey(0.8, 0.8))

		# }

	dfhabs <- data.frame(lv_input$y/375) %>%
		tidyr::gather(hab, cover) %>%
		right_join(dplyr::distinct(dplyr::select(dffilt, hab, cat2)))
	dfhabs$distkm <- ((lv_input$mindists * lv_input$sddist) + lv_input$mndist)/1000

	ggplot(dffilt, aes(x = distkm, y = X50.)) +
	geom_point(data = dfhabs, aes(y = cover), size = 0.5, color = "grey20", shape = 1) +
	geom_line() +
	facet_grid( ~ cat2)+
	geom_ribbon(aes(ymin = X2.5., ymax = X97.5., linetype = NA), colour = 'grey80', alpha = 0.3) +
	theme_classic() +
	xlab('Distance to nearest log ponds (km)') +
	ylab('Cover (proportion)')

}
