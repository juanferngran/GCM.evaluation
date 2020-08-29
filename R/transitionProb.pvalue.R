#' @title Proportion Z-test for WT transition probabilities
#' @description A wrapper of prop.test to assess the significance of WT transition probability between GCM and reanalysis
#' @param obs.grid Cluster grid of observations (reanalysis)
#' @param gcm.grid Cluster grid of GCM simulation 
#' @param what Parameter to be returned. Default to \code{"p.value"}, returning the p-value of the test.
#' See \code{\link[stats]{prop.test}} return values, for the meaning of the other possible values, namely
#' \code{"statistic"}, \code{"parameter"} and \code{"estimate"}.  
#' @param ... Further arguments passed to \code{\link[stats]{prop.test}}
#' @details The function receives as input cluster grids, as returned by \code{\link[transformeR]{clusterGrid}}.
#' 
#' The null hypothesis H0 is that two populations from which the WT transitions were drawn have the same true
#'  proportion of such transition (i.e., for p-value > 0.05 at a 95% confidence interval, H0 can't be rejected).
#' The alternative hypothesis A is that this proportion is different in at least one of the populations
#' (i.e., the GCM has a different distribution).
#' 
#' @return a numeric square matrix with the proportion Z-test results for each transition probability
#' 
#' @references How the test works, it is well explained here: \url{http://www.sthda.com/english/wiki/two-proportions-z-test-in-r}
#' @importFrom transformeR getWT
#' @importFrom stats prop.test
#' @importFrom magrittr %>% extract2
#' 
#' @examples 
#' load("./Data/WTs_cmip5.RData", verbose = TRUE)
#' load("./Data/clustering_eraInterim.RData", verbose = TRUE)

#' pvalues <- transitionProb.test(obs.grid = clusters.era.interim,
#'                                gcm.grid = wts.EC.earth)
#' 
#' pcolors <- RColorBrewer::brewer.pal(n = 9, "YlOrRd") %>% colorRampPalette()
#' # NOTE that levelplot will depict the transposed transition matrix by default!
#' 
#' 
#' lattice::levelplot(t(pvalues), col.regions = pcolors(201), at = seq(0,1,0.01),
#'                    main = "EC-Earth prop.test - p.values",
#'                    xlab = "", ylab = "",
#'                    scales = list(alternating = 3))
#' 

transitionProb.test <- function(obs.grid, gcm.grid, what = "p.value", ...) {
    arg.list  <- list(...)
    wts.obs <- getWT(obs.grid)
    wts.gcm <- getWT(gcm.grid)
    if (is.null(wts.obs) | is.null(wts.gcm)) stop("Input is not a cluster grid")
    obs.grid <- gcm.grid <- NULL
    what <- match.arg(what, choices = c("statistic", "parameter", "p.value", "estimate"))
    wto1 <- wts.obs[-length(wts.obs)]
    wto2 <- wts.obs[-1]
    wtp1 <- wts.gcm[-length(wts.gcm)]
    wtp2 <- wts.gcm[-1]
    # compute all transition probabilities as a vector (cross.freqs*)
    tlabels <- paste(names(wto1), names(wto2), sep = "-->")
    cross.freqs.obs <- table(tlabels) 
    wtnames.obs <- unique(names(wts.obs))
    
    tlabels <- paste(names(wtp1), names(wtp2), sep = "-->")
    cross.freqs.gcm <- table(tlabels) 
    wtnames.gcm <- unique(names(wts.gcm))
    
    # Fixed WT ordering (in decreasing order of annual frequency)
    wtnames.ref <- c("A", "C", "SW", "W", "AW", "NW", "S", "ASW", "N", "ANW",
                     "SE", "CSW", "CS", "E", "NE", "AS", "CSE", "AN", "CW",
                     "CN", "CNW", "ASE", "CNE", "CE", "AE", "ANE")
    ind.order <- match(wtnames.ref, wtnames.obs) %>% na.omit()
    wtnames.obs <- wtnames.obs[ind.order]
    
    ind.order <- match(wtnames.ref, wtnames.gcm) %>% na.omit()
    wtnames.gcm <- wtnames.gcm[ind.order]
    
    # Include missing transitions with zero probability in observations
    ref <- expand.grid(wtnames.ref, wtnames.ref)
    allcombs <- paste(ref[ ,1], ref[ ,2], sep = "-->")
    not.in.obs <- which(!(allcombs %in% names(cross.freqs.obs)))
    namesaux.obs <- allcombs[not.in.obs]
    aux <- rep(0, length(namesaux.obs))
    names(aux) <- namesaux.obs
    all.freqs.obs <- append(cross.freqs.obs, aux)
    
    # Include missing transitions with zero probability in GCM
    # ref <- expand.grid(wtnames.gcm, wtnames.gcm)
    # allcombs <- paste(ref[ ,1], ref[ ,2], sep = "-->")
    not.in.gcm <- which(!(allcombs %in% names(cross.freqs.gcm)))
    namesaux.gcm <- allcombs[not.in.gcm]
    aux <- rep(0, length(namesaux.gcm))
    names(aux) <- namesaux.gcm
    all.freqs.gcm <- append(cross.freqs.gcm, aux)
    
    # Compute matrix of proportion Z-test results for all possible transitions
    mat <- sapply(1:length(wtnames.ref), function(i) {
        wt <- wtnames.ref[i]
        all.freqs.from.obs <- all.freqs.obs[grep(paste0("^", wt, "-->"), names(all.freqs.obs))]
        all.freqs.from.gcm <- all.freqs.gcm[grep(paste0("^", wt, "-->"), names(all.freqs.gcm))]
        nobs <- sum(all.freqs.from.obs)
        nsim <- sum(all.freqs.from.gcm)
        arg.list[["n"]] <- c(nobs, nsim)
        a <- sapply(names(all.freqs.from.obs), USE.NAMES = TRUE, function(j) {
            # tofreq <- ifelse(is.na(all.freqs.from.gcm[j]), 0, all.freqs.from.gcm[j])
            arg.list[["x"]] <- c(all.freqs.from.obs[j], all.freqs.from.gcm[j])
            tryCatch(expr = {
                suppressWarnings({# NaN is returned when the 'to' state does not occur
                    do.call("prop.test", args = arg.list) %>% extract2(what)
                })
            }, error = function(err) NaN)
        })
        ind <- match(paste(wt, wtnames.ref, sep = "-->"), names(a))
        # print(ind)
        return(a[ind])
    })
    mat <- t(mat)
    colnames(mat) <- rownames(mat) <- wtnames.ref
    return(mat)
}

