# #Divergencias de Kullback-Leibler 
# #CMIP5
# erain.WTs <- freqWT(clusters.era.interim)/100
# grid.list <- list(wts.cccma, wts.cnrm, wts.EC.earth, wts.noaa.gfdl,
#                   wts.mohc, wts.ipsl, wts.miroc, wts.mpi.esm.lr, wts.ncc.nor)
# l <- length(grid.list)
# 
# #wts index
# wt.index <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, CMIP5.names))
# for (i in 1:l){
#   wt.index[ ,i] <- freqWT(grid.list[[i]])
# }
# freqs.WTs <- wt.index/100
# 
# divergence <- vector("numeric", length = l)
# for (i in 1:l){
#   divergence[i] <- sum(freqs.WTs[,i]*log(freqs.WTs[,i]/erain.WTs))
# }
# names(divergence) <- CMIP5.names
# sort(divergence)
# 
# #CMIP6
# grid.list <- list(wts.cccma.cmip6, wts.cnrm.cmip6, wts.ec_earth.cmip6.cmip6, wts.gfdl.cmip6,
#                   wts.ukesm1.cmip6, wts.ipsl.cmip6, wts.miroc.cmip6, wts.mpi.lr.cmip6, wts.ncc.nor.cmip6)
# l <- length(grid.list)
# 
# #wts index
# wt.index <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, CMIP6.names))
# for (i in 1:l){
#   wt.index[ ,i] <- freqWT(grid.list[[i]])
# }
# freqs.WTs.cmip6 <- wt.index/100
# 
# divergence.2 <- vector("numeric", length = l)
# for (i in 1:l){
#   divergence.2[i] <- sum(freqs.WTs.cmip6[,i]*log(freqs.WTs.cmip6[,i]/erain.WTs))
# }
# names(divergence.2) <- names.cmip6
# sort(divergence.2)

#Test Divergencias de Kullback-Leibler con funciones de R:

reanalysis.names <- c("ERA-Interim","ERA-20C", "NCEP")

  erain.WTs <- getLWT(clusters.jra) 
  erain.WTs.DJF <- subsetGrid(clusters.jra, season = DJF) %>% getLWT() 
  erain.WTs.MAM <- subsetGrid(clusters.jra, season = MAM) %>% getLWT() 
  erain.WTs.JJA <- subsetGrid(clusters.jra, season = JJA) %>% getLWT() 
  erain.WTs.SON <- subsetGrid(clusters.jra, season = SON) %>% getLWT() 
  
  grid.list.reanalysis <- list(clusters.era.interim, clusters.era20, clusters.ncep)
  grid.list.cmip5 <- list(wts.cccma, wts.cnrm, wts.EC.earth, wts.noaa.gfdl,
                    wts.mohc, wts.ipsl, wts.miroc, wts.mpi.esm.lr, wts.ncc.nor)
  grid.list.cmip6 <- list(wts.cccma.cmip6, wts.cnrm.cmip6, wts.ec_earth.cmip6.cmip6, wts.gfdl.cmip6,
                    wts.ukesm1.cmip6, wts.ipsl.cmip6, wts.miroc.cmip6, wts.mpi.lr.cmip6, wts.ncc.nor.cmip6)
  n <- length(grid.list.reanalysis)
  l <- length(grid.list.cmip5)
  
  #GET LWTs:
  #Reanalysis:
  wts.reanalysis <- matrix(nrow = 26, ncol = n, byrow = TRUE, dimnames = list(wt.names, reanalysis.names))
  wts.reanalysis.DJF <- matrix(nrow = 26, ncol = n, byrow = TRUE, dimnames = list(wt.names, reanalysis.names))
  wts.reanalysis.MAM <- matrix(nrow = 26, ncol = n, byrow = TRUE, dimnames = list(wt.names, reanalysis.names))
  wts.reanalysis.JJA <- matrix(nrow = 26, ncol = n, byrow = TRUE, dimnames = list(wt.names, reanalysis.names))
  wts.reanalysis.SON <- matrix(nrow = 26, ncol = n, byrow = TRUE, dimnames = list(wt.names, reanalysis.names))
  
  for (i in 1:n){
    wts.reanalysis[ ,i] <- getLWT(grid.list.reanalysis[[i]]) 
    
    wts.reanalysis.DJF[ ,i] <- subsetGrid(grid.list.reanalysis[[i]], season = DJF) %>% getLWT() 
    wts.reanalysis.MAM[ ,i] <- subsetGrid(grid.list.reanalysis[[i]], season = MAM) %>% getLWT() 
    wts.reanalysis.JJA[ ,i] <- subsetGrid(grid.list.reanalysis[[i]], season = JJA) %>% getLWT() 
    wts.reanalysis.SON[ ,i] <- subsetGrid(grid.list.reanalysis[[i]], season = SON) %>% getLWT() 
    
  }
  
  #CMIP5:
  wts.cmip5 <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, CMIP5.names))
  wts.cmip5.DJF <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, CMIP5.names))
  wts.cmip5.MAM <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, CMIP5.names))
  wts.cmip5.JJA <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, CMIP5.names))
  wts.cmip5.SON <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, CMIP5.names))
  
  for (i in 1:l){
    wts.cmip5[ ,i] <- getLWT(grid.list.cmip5[[i]]) 
    
    wts.cmip5.DJF[ ,i] <- subsetGrid(grid.list.cmip5[[i]], season = DJF) %>% getLWT() 
    wts.cmip5.MAM[ ,i] <- subsetGrid(grid.list.cmip5[[i]], season = MAM) %>% getLWT() 
    wts.cmip5.JJA[ ,i] <- subsetGrid(grid.list.cmip5[[i]], season = JJA) %>% getLWT() 
    wts.cmip5.SON[ ,i] <- subsetGrid(grid.list.cmip5[[i]], season = SON) %>% getLWT() 
    
  }

  #CMIP6:
  wts.cmip6 <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, names.cmip6))
  wts.cmip6.DJF <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, names.cmip6))
  wts.cmip6.MAM <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, names.cmip6))
  wts.cmip6.JJA <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, names.cmip6))
  wts.cmip6.SON <- matrix(nrow = 26, ncol = l, byrow = TRUE, dimnames = list(wt.names, names.cmip6))
  
  for (i in 1:l){
    wts.cmip6[ ,i] <- getLWT(grid.list.cmip6[[i]])
    
    wts.cmip6.DJF[ ,i] <- subsetGrid(grid.list.cmip6[[i]], season = DJF) %>% getLWT() 
    wts.cmip6.MAM[ ,i] <- subsetGrid(grid.list.cmip6[[i]], season = MAM) %>% getLWT() 
    wts.cmip6.JJA[ ,i] <- subsetGrid(grid.list.cmip6[[i]], season = JJA) %>% getLWT() 
    wts.cmip6.SON[ ,i] <- subsetGrid(grid.list.cmip6[[i]], season = SON) %>% getLWT() 
    
  }
  
  #Divergence:
  
  #DJF: 
  vector1 <- cbind(wts.reanalysis.DJF, wts.cmip5.DJF,wts.cmip6.DJF) %>% t()
  vector2 <- t(erain.WTs.DJF)
  x <- rbind(vector1, vector2)
  
  divergence.DJF <- philentropy::KL(x, unit = "log",  est.prob = "empirical")[nrow(x),1:(nrow(x)-1)]
  names(divergence.DJF) <- c(reanalysis.names, CMIP5.names, names.cmip6)
  #sort(divergence.DJF)
  
  #MAM: 
  vector1 <- cbind(wts.reanalysis.MAM, wts.cmip5.MAM,wts.cmip6.MAM) %>% t()
  vector2 <- t(erain.WTs.MAM)
  x <- rbind(vector1, vector2)
  
  divergence.MAM <- KL(x, unit = "log",  est.prob = "empirical")[nrow(x),1:(nrow(x)-1)]
  names(divergence.MAM) <- c(reanalysis.names, CMIP5.names, names.cmip6)
  
  #JJA: 
  vector1 <- cbind(wts.reanalysis.JJA, wts.cmip5.JJA,wts.cmip6.JJA) %>% t()
  vector2 <- t(erain.WTs.JJA)
  x <- rbind(vector1, vector2)
  
  divergence.JJA <- KL(x, unit = "log",  est.prob = "empirical")[nrow(x),1:(nrow(x)-1)]
  names(divergence.JJA) <- c(reanalysis.names, CMIP5.names, names.cmip6)
  
  #SON: 
  vector1 <- cbind(wts.reanalysis.SON, wts.cmip5.SON,wts.cmip6.SON) %>% t()
  vector2 <- t(erain.WTs.SON)
  x <- rbind(vector1, vector2)
  
  divergence.SON <- KL(x, unit = "log",  est.prob = "empirical")[nrow(x),1:(nrow(x)-1)]
  names(divergence.SON) <- c(reanalysis.names, CMIP5.names, names.cmip6)
  
  #Year: 
  vector1 <- cbind(wts.reanalysis, wts.cmip5, wts.cmip6) %>% t()
  vector2 <- t(erain.WTs)
  x <- rbind(vector1, vector2)
  
  divergence.Year <- KL(x, unit = "log",  est.prob = "empirical")[nrow(x),1:(nrow(x)-1)]
  names(divergence.Year) <- c(reanalysis.names, CMIP5.names, names.cmip6)

  #Create KLdivergence.DF:
  KLdivergence.DF <- data.frame(GCM = factor(rep(c(reanalysis.names, CMIP5.names, CMIP6.names),5), 
                                             levels = c(reanalysis.names, CMIP5.names[10:1])),
                            Experiment = factor(rep(c(rep("Reanalysis", 3), rep("CMIP5", 9), rep("CMIP6", 9)),5),
                                                levels = c("Reanalysis", "CMIP5","CMIP6")),
                            Season = factor(c(rep("DJF",21),rep("MAM",21),rep("JJA",21),rep("SON",21), rep("Year",21)), 
                                            levels = order.season),
                            KL = c(as.vector(divergence.DJF), as.vector(divergence.MAM), as.vector(divergence.JJA),
                                     as.vector(divergence.SON), as.vector(divergence.Year)))
  KLdivergence.DF <- KLdivergence.DF[order(KLdivergence.DF$Experiment),]
  row.names(KLdivergence.DF) <- 1:105
  
  save(KLdivergence.DF, file = "KLdivergence.DF.RData")
  write.csv(KLdivergence.DF, file = "KLdivergence.DF.csv")

  #PLOT:
  
  dev.new()
  pcolors <- RColorBrewer::brewer.pal(n = 9, "OrRd") %>% colorRampPalette()
  
  
  z <- subset(KLdivergence.DF, Experiment == "Reanalysis")
  p1 <- lattice::levelplot(KL ~ Season + GCM , data = z, 
                           col.regions = pcolors(201), set.min=0, set.max=0.2,
                           at = seq(0, .22, .01),
                           colorkey = list(space = 'bottom'),
                           par.settings=list(layout.widths=list(key.ylab.padding = 2))) +  
    lattice::xyplot(KL ~ Season + GCM , data = z,
                    panel = function(x, y, ...) {
                      ltext(x = z$Season, y = z$GCM, labels = round(z$KL, 3), cex = 0.9, font = 1,
                            fontfamily = "HersheySans")})
  
  w <- subset(KLdivergence.DF, Experiment == "CMIP5")
  p2 <- lattice::levelplot(KL ~ Season + GCM , data = w, 
                           col.regions = pcolors(201), set.min=0, set.max=0.2,
                           at = seq(0, .22, .01),
                           scales=list(alternating=1),
                           colorkey = list(space = 'bottom')) +  
    lattice::xyplot(KL ~ Season + GCM , data = w,
                    panel = function(x, y, ...) {
                      ltext(x = w$Season, y = factor(CMIP5.names, levels = CMIP5.names[9:1]), labels = round(w$KL, 3), cex = 0.9, font = 1,
                            fontfamily = "HersheySans")})
  
  v <- subset(KLdivergence.DF, Experiment == "CMIP6")
  p3 <- lattice::levelplot(KL ~ Season + GCM , data = v, 
                           col.regions = pcolors(201), set.min=0, set.max=0.2,
                           at = seq(0, .22, .01),
                           scales=list(y = list(alternating=2)),
                           colorkey = list(space = 'bottom')) +  
    lattice::xyplot(KL ~ Season + GCM , data = v,
                    panel = function(x, y, ...) {
                      ltext(x = v$Season, y = factor(CMIP5.names, levels = CMIP5.names[9:1]), labels = round(v$KL, 3), cex = 0.9, font = 1,
                            fontfamily = "HersheySans")})
  
  comb_levObj <- c(Reanalysis = p1, CMIP5 = p2, CMIP6 = p3, layout = c(3, 1))
  print(comb_levObj)
  update(comb_levObj, main = "Kullback-Leibler Divergence: GCM & Reanalysis vs NCEP", xlab = NULL,
         scales = list(x = list(alternating=3), y = list(rot=0)))
  
  
  
              