#' Create flow diagram of model
#'
#' @author Christopher J. Brown
#' @rdname habitat_loadings_plot
#' @export


habitat_loadings_plot <- function(file, prows=2, pcols=2){
data(lv_input)
load(file = file)

smcmn <- summary(mcout3)


#Weights on habitats

dimnames(lv_input$y)[[2]]
habnams <- dimnames(lv_input$y)[[2]]
habnams_full <- lv_input$habnams_full$CATEGORY
nhabs <- length(habnams)

#Number of LVs
colnams <- dimnames(mcout3[[1]])[[2]]
ialpha <- stringr::str_detect(colnams, 'alpha')
num_levels <- sum(ialpha)/nhabs

par(mfrow = c(prows, pcols), mar = c(5,6,2,2))
iord <- plot_beta_mnm(smcmn, habnams_full, xlab = "Constrained LV")
par(mar = c(5,0.5,2,0.5))
purrr::walk(1:num_levels, ~plot_alpha_mnm(smcmn, habnams_full, lv =.x, nlv = num_levels, iord, xlab = paste("LV", .x)))



}
