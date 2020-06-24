
### Helpers

getLWT <- function(grid){
  wt.names <-c("A", "ANE", "AE", "ASE", "AS", "ASW", "AW", "ANW", "AN",
               "NE",  "E", "SE",  "S",  "SW",  "W",  "NW",  "N", 
               "C", "CNE", "CE", "CSE", "CS", "CSW", "CW", "CNW", "CN")
  
  out <- attr(x = grid, which = "wt.index") %>% table()
  out <- out[match(1:26, names(out))]
  out[which(is.na(out))] <- 0
  names(out) <- wt.names
  return(out)
} 

###########

loadCMIP5 <- function(historical = NULL,
                      rcp8.5 = NULL,
                      var = "psl",
                      season = NULL,
                      lonLim = c(-45, 66),
                      latLim = c(22,73)){
  
  grid1 <- loadGridData(dataset = historical,
                        var = var,
                        lonLim = lonLim,
                        latLim = latLim,
                        years = 1980:2005,
                        time = "DD", 
                        aggr.d = "mean")
  
  grid2 <- loadGridData(dataset = rcp8.5,
                        var = var,
                        lonLim = lonLim,
                        latLim = latLim,
                        years = 2006:2010,
                        time = "DD", 
                        aggr.d = "mean")
  out <- bindGrid(grid1, grid2, dimension = "time")
  out <- subsetGrid(out, season = season)
  if(getShape(out, dimension = "time") != 10957){
    warning("Output grid has not 10957 days.")
  }
  return(out)
}

#########

freqWT <- function(grid = NULL, 
                   season = NULL){
  
  if(!is.null(season)){
    grid <- subsetGrid(grid, season = season)
  }
  wt.names <-c("A", "ANE", "AE", "ASE", "AS", "ASW", "AW", "ANW", "AN",
                     "NE",  "E", "SE",  "S",  "SW",  "W",  "NW",  "N", 
               "C", "CNE", "CE", "CSE", "CS", "CSW", "CW", "CNW", "CN")
  
  out <- getWT(grid) %>% names() %>% table() %>% prop.table() 
  out <- out[match(wt.names, names(out))]
  out[which(is.na(out))] <- 0
  names(out) <- wt.names
  out <- round(out*100, 2)
  return(out)
}

#########

relativeFreqWT <- function(..., ref, 
                           cluster = 1:26,
                           season = NULL){
  #browser()
  grid.list <- list(...)
  l <- length(grid.list)
  
  #wts index
  wt.index <- matrix(nrow = 26, ncol = l, byrow = TRUE)
  for (i in 1:l){
    wt.index[ ,i] <- freqWT(grid.list[[i]], season = season)
  }
  
  ### Diff. among GCM and Reanalysis: 
  
  diff.freqs <- matrix(nrow = length(cluster), ncol = l, byrow = TRUE)
  for (i in 1:l){
    diff.freqs[ ,i] <- wt.index[cluster,i] - ref[cluster]
  }
  return(diff.freqs)
}

#########

getLambWTIndex <- function(){
  
  
}

mse <- function(data) { mean(data^2) }

rmse <- function(data){ sqrt(mean(data^2)) }

divergence <- function(data, ref){ sum(freqs.WTs[,i]*log(freqs.WTs[,i]/erain.WTs)) }


freqWT.pvalues <- function(data.mat = NULL, #Matriz donde se colocarán las cruzes
                           season = NULL,
                           clusters = NULL, 
                           ref = NULL, 
                           alpha = 0.05){
  #     browser()
  
  if(!is.null(season)){
    ref <- subsetGrid(ref, season = season)
  }
  
  names.cmip5 <- c("CanESM2", "CNRM-CM5", "EC-EARTH", "GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5-LR", "MIROC5", "MPI-ESM-LR", "NorESM1-M")
  names.cmip6 <- c("CanESM5", "CNRM-CM6-1","EC-EARTH3", "GFDL-ESM4", "UKESM1-0-LL", "IPSL-CM6A-LR", "MIROC6",  "MPI-ESM1-2-LR", "NorESM2-LM")
  GCM.names <- c("JRA", "ERA-20C", "NCEP", names.cmip5, names.cmip6)
  models <- c("clusters.jra","clusters.era20","clusters.ncep",
              "wts.cccma", "wts.cnrm", "wts.EC.earth", "wts.noaa.gfdl",
              "wts.mohc", "wts.ipsl", "wts.miroc", "wts.mpi.esm.lr","wts.ncc.nor",
              "wts.cccma.cmip6", "wts.cnrm.cmip6", "wts.ec_earth.cmip6.cmip6",
              "wts.gfdl.cmip6", "wts.ukesm1.cmip6", "wts.ipsl.cmip6",
              "wts.miroc.cmip6", "wts.mpi.lr.cmip6", "wts.ncc.nor.cmip6")
  
  season.freqs.matrix <- matrix(nrow = length(models), ncol = length(clusters), 
                                dimnames = list(GCM.names, clusters))
  
  wt.names <-c("A", "ANE", "AE", "ASE", "AS", "ASW", "AW", "ANW", "AN",
               "NE",  "E", "SE",  "S",  "SW",  "W",  "NW",  "N", 
               "C", "CNE", "CE", "CSE", "CS", "CSW", "CW", "CNW", "CN")
  WTS <- match(clusters, wt.names)
  
  for (i in 1:length(models)) {
    model <- get(models[i])
    if(!is.null(season)){
      model <- subsetGrid(model, season = season)
      season.freqs.matrix[i,] <- table(getWT(model))[WTS]
    }else{
      season.freqs.matrix[i,] <- table(getWT(model))[WTS]
    }
  }
  
  arg.list <- list()
  mat <- sapply(1:dim(season.freqs.matrix)[1], function(i) {
    freqs.obs <- table(getWT(ref))[WTS] 
    freqs.sim <- season.freqs.matrix[i, ]
    nobs <- sum(freqs.obs)
    nsim <- sum(freqs.sim)
    arg.list[["n"]] <- c(nobs, nsim)
    z <- sapply(1:length(freqs.obs), function(j) {
      arg.list[["x"]] <- c(freqs.obs[[j]], freqs.sim[[j]])
      suppressWarnings({
        do.call("prop.test", args = arg.list) %>% extract2("p.value")
      })
    })
    return(z)
  })
  
  ind <- match(colnames(data.mat), rownames(season.freqs.matrix))
  
  # Significantly different transition probabilities
  # Filter those whose frequency is zero either in the model or in the reanalysis 
  # (these are always signif)
  signif <- which((mat[ ,ind] < alpha), arr.ind = TRUE)
  points <- matrix(signif, nrow(signif))
  return(points)
  
}  
