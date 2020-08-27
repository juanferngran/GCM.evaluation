#' @title Transition probability matrix score (TPMS)
#' @description Provide an overall measure summarizing transition probability model performance. See details.
#' @param obs.grid Cluster grid of observations (reanalysis)
#' @param gcm.grid Cluster grid of GCM simulation 
#' @param conf.level confidence level. Default to 0.95
#' @param include.nonexisting Logical flag. Should missing transitions be included to compute the score?
#'  Default to \code{FALSE}. See Details.
#' @return A numeric value for the calculated TPMS (0 is perfect)
#' @importFrom magrittr %>% 
#' @export
#' @author juaco
#' @details The TPMS It is intended to serve as a summary measure that allows to rank different GCMs according
#'  to their transition probability matrix "fingerprint".
#'  
#' By default (\code{include.nonexisting = FALSE}), the summary measure is the summation of the absolute differences
#' between model and reanalysis transition probabilities, considering only those model transitions that are
#' significantly different from the observation (reanalysis). This choice will only consider the transition
#' probabilities that exist bot in the GCM and the reanalysis, thus not accounting for the "missing" transitions, 
#' i.e.: either transitions that exist in the reanalysis but are never simulated by the model,
#' or transitions that are simulated by the model but never happen in the reanalysis.
#' 
#' In order to account for these missing transitions, the argument \code{include.nonexisting} is set to \code{TRUE}.
#' In this case, the missing probability (either from GCM or for reanalysis) is assigned a value of zero, and included
#' in the TMPS calculation. Thus, this version of TMPS will be always equal or greater than the default one.
#' 
#' @examples \dontrun{
#' load("./juan_WORK/Data/WTs_cmip5_europe.RData", verbose = TRUE)
#' TPMS(wts.EC.earth, clusters.era.interim)
#' # This is the preliminary choice for the Lamb paper
#' TPMS(wts.EC.earth, clusters.era.interim, include.nonexisting = TRUE)
#' }

# source("R/transitionProb.R")
# source("R/transitionProb.pvalue.R")

TPMS <- function(obs.wt.grid,
                 gcm.wt.grid,
                 conf.level = .05,
                 include.nonexisting = TRUE) {
    stopifnot(is.logical(include.nonexisting))
    tprob.matrix = transitionProb(gcm.wt.grid)
    pval.matrix = transitionProb.test(obs.grid = obs.wt.grid, gcm.grid =  gcm.wt.grid, conf.level = conf.level)
    tprob.ref = transitionProb(obs.wt.grid)
    gcm.wt.grid <- obs.wt.grid <- NULL
    sig <- new <- missing <- 0 
    # Significantly different transition probabilities
    signif <- which((pval.matrix < conf.level & !is.na(tprob.matrix) & !is.na(tprob.ref)), arr.ind = TRUE)
    if (nrow(signif) > 0L) {
        sig <- sapply(1:nrow(signif), FUN = function(i) {
            p_gcm <- tprob.matrix[signif[i,1], signif[i,2]]
            p_rea <- tprob.ref[signif[i,1], signif[i,2]]
            abs(p_gcm - p_rea)
        }) %>% sum()
    }
    if (include.nonexisting) {
        # Transitions that do not exist in reanalysis, but exist in the model 
        newt <- which((is.na(tprob.ref) & !is.na(tprob.matrix)), arr.ind = TRUE)
        if (nrow(newt) > 0L) {
            new <- sapply(1:nrow(newt), FUN = function(i) {
                tprob.matrix[newt[i,1], newt[i,2]]
            }) %>% sum()
        }
        # Transitions that exist in reanalysis, but not in the model 
        missingt <- which((!is.na(tprob.ref) & is.na(tprob.matrix)), arr.ind = TRUE)
        if (nrow(missingt) > 0L) {
            missing <- sapply(1:nrow(missingt), FUN = function(i) {
                tprob.ref[missingt[i,1], missingt[i,2]]
            }) %>% sum()
        } 
    }
    sum(sig, new, missing) %>% return()
}


