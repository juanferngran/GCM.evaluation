#' @title Transition probability model assessment
#' @description A wraper summarizing transition probability results from GCMs,
#'  as compared to a given reference.
#' @param tprob.matrix GCM transition probability matrix, as returned
#'   by \code{transitionProb}
#' @param pval.matrix Optional. p-value matrix, as returned by \code{transitionProb.pvalue}
#' @param tprob.ref Optional (required if \code{pval.matrix} is not \code{NULL}).
#'  Transition probability object from observed reference, 
#'  as returned by \code{transitionProb}
#' @param alpha alpha cut-off value. Default to 0.05 (i.e. 95% confidence interval)
#' @param title Graph title
#' @details If \code{pval.matrix} is set to \code{NULL}, 
#' @return A transition probability asessment plot
#' @importFrom graphics image legend axis
#' @importFrom fields image.plot

tprobPlot <- function(tprob.matrix,
                      pval.matrix = NULL,
                      tprob.ref = NULL,
                      alpha = 0.05,
                      title = NULL) {
    # if (is.null(pcolors)) {
    pcolors <- RColorBrewer::brewer.pal(n = 9, "YlOrRd")[-1] %>% colorRampPalette()
    # }
    wtnames <- colnames(tprob.matrix)
    par(mar = c(9, 4.5, 9, 9))
    image(t(tprob.matrix), col = pcolors(201), axes = FALSE, main = title)
    axis(1, at = seq(0, 1, 1/25), labels = wtnames, las = 2)
    axis(2, at = seq(0, 1, 1/25), labels = wtnames, las = 1)
    axis(3, at = seq(0, 1, 1/25), labels = wtnames, las = 2)
    axis(4, at = seq(0, 1, 1/25), labels = wtnames, las = 1)
    fields::image.plot(tprob.matrix, legend.only = TRUE, col = pcolors(21))

    if (!is.null(pval.matrix)) {
        stopifnot(!is.null(tprob.ref))
        
        coords <- seq(0, 1, 1/25)
        
        # Empty transitions (for testing purposes)
        # empty.model <- which(is.na(tprob.matrix), arr.ind = TRUE)
        # for (i in 1:nrow(empty.model)) {
        #     points(coords[empty.model[i,2]], coords[empty.model[i,1]], pch = "/")
        # }
        
        # Transitions that do not exist in reality, but exist in the model 
        newt <- which((is.na(tprob.ref) & !is.na(tprob.matrix)), arr.ind = TRUE)
        for (i in 1:nrow(newt)) {
            points(coords[newt[i,2]], coords[newt[i,1]], pch = 21)
        }
        
        # Transitions that do not exist in model, but exist in reality 
        missingt <- which((!is.na(tprob.ref) & is.na(tprob.matrix)), arr.ind = TRUE)
        for (i in 1:nrow(missingt)) {
            points(coords[missingt[i,2]], coords[missingt[i,1]], pch = 19, cex = .5)
        }
        
        # Significantly different transition probabilities
        # Filter those whose frequency is zero either in the model or in the reanalysis 
        # (these are always signif)
        signif <- which((pval.matrix < alpha & !is.na(tprob.matrix) & !is.na(tprob.ref)), arr.ind = TRUE)
        # signif <- which(pval.matrix < alpha, arr.ind = TRUE)
        
        if (nrow(signif) > 0L) {
            for (i in 1:nrow(signif)) {
                points(coords[signif[i,2]], coords[signif[i,1]], pch = 4, cex = 1.5)
            }
        }
        
        # Legend
        
        par(xpd = TRUE) # Allow putting legend outside graphic area
        legend(0.1, -0.15,
               c("Does not occur in reanalysis", 
                 "Does not occur in the model (but it does in reanalysis)",
                 "Significantly different probability"),
               pch = c(21, 19, 4), pt.cex = c(1, 0.5, 1.5),
               cex = .85, bty = "n")
    }
}


    
