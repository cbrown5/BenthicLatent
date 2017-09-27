#' Extract sample parameters from an LV model
#' 
#' Only extracts parameters for a single chain 
#' @author Christopher J. Brown
#' @rdname extract_params
#' @export


extract_params <- function(mcout){
	
	colnams <- dimnames(mcout[[1]])[[2]]

	ialpha <- str_detect(colnams, 'alpha')
	ibeta <- str_detect(colnams, 'beta')
	iint <- str_detect(colnams, 'int')
	igam <- str_detect(colnams, '\\bgamma\\b')
	igamf <- str_detect(colnams, 'gammaf')
	inu <- str_detect(colnams, '\\bnu\\b')
	inuwq <- str_detect(colnams, '\\bnuwq\\b')

	mcmat <- as.matrix(mcout[[1]])

	params <- list(alpha = mcmat[,ialpha], 
	beta = mcmat[,ibeta], 
	int = mcmat[,iint],
	gamma = mcmat[,igam],
	gammaf = mcmat[,igamf], 
	nu = mcmat[,inu], 
	nuwq = mcmat[,inuwq], 
	nsamp = nrow(mcmat))

	return(params)
	}
