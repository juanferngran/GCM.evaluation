#' @title Transition probability matrix
#' @description Compute a Lamb WT transition probability matrix
#' @param grid A grid with a (named) \code{wt.index} attribute (as returned by \code{\link[transformeR]{lambWT}})
#' @param missing.to.zero When a certain transition does never occur, should it be set to zero probability, or indicated as NA?. Default to FALSE (i.e., NAs are preserved).
#' @return A 2D matrix of transition probabilities. The matrix must be read in the form "row-->column" as the transition probability "from-->to"
#' @details The frequency of transition is calculated as the transition probability from WT-A to WT-B
#' divided by the frequency of all transitions from WT-A to all other WTs
#' @importFrom transformeR getWT
#' @author juaco
#' @examples \dontrun{
#' load("/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/clustering_eraInterim.RData", verbose = TRUE)
#' 
#' Transition probability matrix ----------------------
#' mat <- transitionProb(clusters.era.interim)
#' # Plotting -----------------------------------------
#' pcolors <- RColorBrewer::brewer.pal(n = 9, "YlOrRd")[-1] %>% colorRampPalette()
#' # NOTE that levelplot will depict the transposed transition matrix by default!
#' lattice::levelplot(t(mat), col.regions = pcolors(201), at = seq(0,.6,0.05),
#'                   main = "Transition probabilities",
#'                   xlab = "", ylab = "",
#'                   scales = list(alternating = 3))
#' }

# Alternative implementation using image.plot (in progress)
# tprob.matrix=mat
# 
# tprobPlot <- function(tprob.matrix, pval.matrix, pcolors = NULL) {
#   if (is.null(pcolors)) {
#     pcolors <- RColorBrewer::brewer.pal(n = 9, "YlOrRd")[-1] %>% colorRampPalette()
#   }
#   wtnames <- colnames(tprob.matrix)
#   par(mar = c(5, 4.5, 4, 9))
#   image(t(tprob.matrix), col = pcolors(201), axes = FALSE) 
#   axis(1, at = seq(0, 1, 1/25), labels = wtnames, las = 2)
#   axis(2, at = seq(0, 1, 1/25), labels = wtnames, las = 1)
#   axis(3, at = seq(0, 1, 1/25), labels = wtnames, las = 2)
#   axis(4, at = seq(0, 1, 1/25), labels = wtnames, las = 1)
#   fields::image.plot(tprob.matrix, legend.only = TRUE, col = pcolors(21))
# }

transitionProb <- function(grid, missing.to.zero = FALSE) {
  wts <- getWT(grid)
  if (is.null(wts)) stop("The \'wt.index\' attribute was not found")
  stopifnot(is.logical(missing.to.zero))
  
  # wt.freqs <- names(wts) %>% table()
  # sort(table(names(wts))/sum(wts), decreasing = T) %>% barplot()
  
  # 1-day lagged series
  wt1 <- wts[-length(wts)]
  wt2 <- wts[-1]
  
  # compute all transition probabilities as a vector
  tlabels <- paste(names(wt1), names(wt2), sep = "-->")
  cross.freqs <- table(tlabels) 
  wtnames <- unique(names(wts))
  
  # Fixed WT ordering (in decreasing order of annual frequency)
  wtnames.ref <- c("A", "C", "SW", "W", "AW", "NW", "S", "ASW", "N", "ANW",
                   "SE", "CSW", "CS", "E", "NE", "AS", "CSE", "AN", "CW",
                   "CN", "CNW", "ASE", "CNE", "CE", "AE", "ANE")
  ind.order <- match(wtnames.ref, wtnames)
  wtnames <- wtnames[ind.order]
  
  # Include missing transitions with zero probability
  ref <- expand.grid(wtnames, wtnames)
  allcombs <- paste(ref[ ,1], ref[ ,2], sep = "-->")
  namesaux <- allcombs[which(!(allcombs %in% names(cross.freqs)))]
  aux <- rep(0, length(namesaux))
  names(aux) <- namesaux
  all.freqs <- append(cross.freqs, aux)
  
  # The frequency of transition is calculated as the transition probability from WT-A to WT-B
  # divided by the frequency of all transitions from WT-A to all other WTs
  
  # Compute matrix of cross-probabilities
  mat <- sapply(1:length(wtnames), USE.NAMES = FALSE, function(i) {
    wt <- wtnames[i]
    freq.from <- all.freqs[grep(paste0("^", wt, "-->"), names(all.freqs))] %>% sum()
    sapply(paste(wt, wtnames, sep = "-->"), function(x) {
      all.freqs[x] / freq.from
    })
  })
  
  # Set non-recorded transitions to missing (to highlight that never occured in the dataset) 
  if (!missing.to.zero) {
    if (isTRUE(any(mat == 0))) mat[which(mat == 0)] <- NA
  }
  mat <- t(mat) 
  colnames(mat) <- rownames(mat) <- wtnames
  return(mat)
}
