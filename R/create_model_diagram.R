#' Create flow diagram of model
#'
#' @author Christopher J. Brown
#' @rdname create_model_diagram
#' @export


create_model_diagram <- function(){

par(mar = c(rep(0.5, 4)))
plot(0,0, xlim = c(0,1), ylim = c(0,1), type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n')

points(0.25, 0.5, cex = 20)
points(0.75, 0.5, cex = 20)

points(0.25, 0.1, cex = 15, pch = 0)
points(0.55, 0.1, cex = 15, pch = 0)

points(0.25, 0.9, cex = 15, pch = 0)
points(0.5, 0.9, cex = 15, pch = 0)
points(0.75, 0.9, cex = 15, pch = 0)

alwd <- 2
alen <- 0.2
acol <- 'grey'

arrows(0.75, 0.65, 0.25, 0.79, lwd = alwd,  len = alen, col = acol)
arrows(0.75, 0.65, 0.55, 0.79, lwd = alwd,  len = alen, col = acol)
arrows(0.75, 0.65, 0.75, 0.79, lwd = alwd,  len = alen, col = acol)

arrows(0.25, 0.65, 0.25, 0.79, lwd = alwd,  len = alen)
arrows(0.25, 0.65, 0.55, 0.79, lwd = alwd,  len = alen)
arrows(0.25, 0.65, 0.75, 0.79, lwd = alwd,  len = alen)

arrows(0.25, 0.21, 0.25, 0.35, lwd = alwd,  len = alen)
arrows(0.55, 0.21, 0.28, 0.35, lwd = alwd,  len = alen)

text(0.25, 0.1, paste('Distance to','\n','nearest sediment','\n','source'), cex = 0.8)
text(0.55, 0.1, 'Flow strength')

text(0.25, 0.5, 'Latent variable 1')
text(0.75, 0.5, 'Latent variable 2')

text(0.25, 0.89, 'Turf algae', font = 3)
text(0.5, 0.89, paste('Branching','\n','coral'), font = 3)
text(0.75, 0.89, 'Silt', font = 3)

}
