# Evaluation of GCMs from CMIP5 and CMIP6:

library(loadeR)
library(transformeR)
library(lattice)
library(latticeExtra)
library(magrittr)

#### Experiment 1 #### 
#Comparación Reanalisis: 

# dataInventory("http://meteo.unican.es/tds5/dodsC/interim/daily/interim20_daily.ncml")

wmo <- 1981:2010
lonLim = c(-45, 66)
latLim = c(22,73)
# jra <- loadGridData(dataset = "http://meteo.unican.es/tds5/dodsC/jra55/daily/JRA55_daily_Dataset.ncml",
#                     var = "slp",
#                     lonLim = lonLim,
#                     latLim = latLim,
#                     season = c(12,1:11),
#                     years = wmo,
#                     time = "DD", 
#                     aggr.d = "mean")
# era.interim <- loadGridData(dataset = "http://meteo.unican.es/tds5/dodsC/interim/daily/interim20_daily.ncml",
#                             var = "SLP",
#                             lonLim = lonLim,
#                             latLim = latLim,
#                             season = c(12,1:11),
#                             years = wmo,
#                             time = "DD", 
#                             aggr.d = "mean")
# 
# 
# #Lamb WTs of JRA and ERA-Interim: 
# clusters.jra <- clusterGrid(jra, type = "lamb")
# clusters.era.interim <- clusterGrid(era.interim, type = "lamb")
# save(era.interim, clusters.era.interim, file = "clustering_eraInterim.RData")
# save(jra, clusters.jra, file = "clustering_jra.RData")

# lonLim = c(-45, 66)
# latLim = c(22,73)

#Lamb Coords (l1):

centerlon = -5
centerlat = 55
lon.array <- rep(centerlon, times=16)+c(-5, 5, -15, -5, 5, 15, -15, -5, 5, 15, -15, -5, 5, 15, -5, 5)
lat.array <- rep(centerlat, times=16)+c(10, 10, 5, 5, 5, 5, 0, 0, 0, 0, -5, -5, -5, -5, -10, -10)
coords <- cbind(lon.array, lat.array) %>% SpatialPoints()
l1 <- list("sp.points", coords, col= 1) #which = 1, pch = 1, lwd = 1.5)

#Load clustering from reanalysis:
load("/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/clustering_eraInterim.RData", verbose = TRUE)
load("/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/clustering_jra.RData", verbose = TRUE)
load("/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/clustering_era20.RData", verbose = TRUE)
load("/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/clustering_ncep.RData", verbose = TRUE)
### Figure 1: BAR PLOT with LWTs comparison of the four reanalysis in all Seasons:
# Subset por estaciones:
DJF <- c(12,1,2)
MAM <- c(3,4,5)
JJA <- c(6,7,8)
SON <- c(9,10,11)
wt.names <-  c("A", "ANE", "AE", "ASE", "AS", "ASW", "AW", "ANW", "AN",
               "NE",  "E", "SE",  "S",  "SW",  "W",  "NW",  "N", 
               "C", "CNE", "CE", "CSE", "CS", "CSW", "CW", "CNW", "CN")

# getWT(clusters.era.interim) %>% names() %>% table() %>% prop.table() 

dev.new()
layout(matrix(c(rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6), 5), ncol=1)) 
par(mai = rep(0.52, 4))

a1 <- freqWT(clusters.era.interim, season = DJF) 
wt.order.a <- sort.int(a1, decreasing = TRUE, index.return = TRUE)$ix
a2 <- freqWT(clusters.jra, season = DJF)[wt.order.a] 
a3 <- freqWT(clusters.ncep, season = DJF)[wt.order.a] 
a4 <- freqWT(clusters.era20, season = DJF)[wt.order.a]
wt.freqs <- matrix(c(a1[wt.order.a],a2,a3,a4), nrow = 4, ncol = 26, byrow = TRUE, 
                   dimnames = list(c("ERA-Interim", "JRA", "NCEP","ERA20"), names(a2)))

p1 <- barplot(wt.freqs, 
              beside = TRUE,  
              col = c("#612C69", "#459ED5", "#FCCF61" , rgb(0.3,0.9,0.4,0.6)), border = "grey",
              ylab = "freq. [%]", ylim = c(0,20),
              main = "DJF")

b1 <- freqWT(clusters.era.interim, season = MAM)
wt.order.b <- sort.int(b1, decreasing = TRUE, index.return = TRUE)$ix
b2 <- freqWT(clusters.jra, season = MAM)[wt.order.b] 
b3 <- freqWT(clusters.ncep, season = MAM)[wt.order.b]  
b4 <- freqWT(clusters.era20, season = MAM)[wt.order.b] 
wt.freqs <- matrix(c(b1[wt.order.b],b2,b3,b4), nrow = 4, ncol = 26, byrow = TRUE, 
                   dimnames = list(c("ERA-Interim", "JRA", "NCEP","ERA20"), names(b2)))
p2 <- barplot(wt.freqs, 
              beside = TRUE,  
              col = c("#612C69", "#459ED5", "#FCCF61" , rgb(0.3,0.9,0.4,0.6)), border = "grey",
              ylab = "freq. [%]", ylim = c(0,20),
              main = "MAM")

c1 <- freqWT(clusters.era.interim, season = JJA) 
wt.order.c <- sort.int(c1, decreasing = TRUE, index.return = TRUE)$ix
c2 <- freqWT(clusters.jra, season = JJA)[wt.order.c] 
c3 <- freqWT(clusters.era20, season = JJA)[wt.order.c] 
c4 <- freqWT(clusters.ncep, season = JJA)[wt.order.c] 
wt.freqs <- matrix(c(c1[wt.order.c],c2,c3,c4), nrow = 4, ncol = 26, byrow = TRUE, 
                   dimnames = list(c("ERA-Interim", "JRA", "NCEP","ERA20"), names(c2)))
p3 <- barplot(wt.freqs, 
              beside = TRUE,  
              col = c("#612C69", "#459ED5", "#FCCF61" , rgb(0.3,0.9,0.4,0.6)), border = "grey",
              ylab = "freq. [%]", ylim = c(0,20),
              main = "JJA")

d1 <- freqWT(clusters.era.interim, season = SON)
wt.order.d <- sort.int(d1, decreasing = TRUE, index.return = TRUE)$ix
d2 <- freqWT(clusters.jra, season = SON)[wt.order.d] 
d3 <- freqWT(clusters.ncep, season = SON)[wt.order.d] 
d4 <- freqWT(clusters.era20, season = SON)[wt.order.d] 
wt.freqs <- matrix(c(d1[wt.order.d],d2,d3,d4), nrow = 4, ncol = 26, byrow = TRUE, 
                   dimnames = list(c("ERA-Interim", "JRA", "NCEP","ERA20"), names(d2)))
p4 <- barplot(wt.freqs, 
              beside = TRUE,  
              col = c("#612C69", "#459ED5", "#FCCF61" , rgb(0.3,0.9,0.4,0.6)), border = "grey",
              ylab = "freq. [%]", ylim = c(0,20),
              main = "SON")

par(mai=c(0,0,0,0))
plot.new()
legend(legend = c("ERA-Interim", "JRA", "NCEP","ERA20") , 
       fill = c("#612C69", "#459ED5", "#FCCF61" , rgb(0.3,0.9,0.4,0.6)), 
       "center", horiz=TRUE, border = "transparent", bty = "n")


### Figure 2: SpatialPlot of annual climatologies from 8-LWTs-subset of ERA-Interim: 

WTs.index <- c(1,18,15,14,16,13,7,17)
subsetWT <- c("A", "C", "W", "SW", "NW", "S", "AW", "N")

t2 <- freqWT(clusters.era.interim)
names.attr2 <- paste0(subsetWT, ": ", t2[WTs.index], "%")

cts2 <- lapply(1:8, function(x) {
  a <- suppressMessages(climatology(subsetGrid(clusters.era.interim, cluster = WTs.index[x])))
})

cts.mg2 <- makeMultiGrid(cts2, skip.temporal.check = TRUE)
dev.new()
breaks <- seq(99300, 103300, 200)
coloykey.labels <- c(99300, 99800, 100300, 100800, 101300, 101800, 102300, 102800, 133000)
visualizeR::spatialPlot(cts.mg2, sp.layout = list(l1), backdrop.theme = "coastline", rev.colors = TRUE,
                        main = "Lamb WTs from ERA-Interim (1981-2010)", useRaster=TRUE,
                        set.min = min(breaks), set.max = max(breaks), at = breaks, 
                        colorkey = list(space = 'bottom',
                                        labels=list(at = seq(99300, 103300, 500), 
                                                    labels = coloykey.labels)),
                        layout = c(2,4), as.table = TRUE, names.attr = names.attr2, contour=TRUE, lty = 3)



####### EXPERIMENT 2 ########

### Load Datasets CMIP5 y merge of 5 years from rcp8.5 as specified in section 2.1 of the study:

# cccma <- loadCMIP5(historical = "http://meteo.unican.es/tds5/dodsC/cmip5/CCCMA/CANESM2/historical/day/cccma_canesm2_historical_r1i1p1.ncml",
#                    rcp8.5 = "http://meteo.unican.es/tds5/dodsC/cmip5/CCCMA/CANESM2/rcp85/day/cccma_canesm2_rcp85_r1i1p1.ncml", season = c(12,1:11))
# cnrm <- loadCMIP5(historical = "http://meteo.unican.es/tds5/dodsC/cmip5/CNRM-CERFACS/CNRM-CM5/historical/day/cnrm-cerfacs_cnrm-cm5_historical_r1i1p1.ncml",
#                   rcp8.5 = "http://meteo.unican.es/tds5/dodsC/cmip5/CNRM-CERFACS/CNRM-CM5/rcp85/day/cnrm-cerfacs_cnrm-cm5_rcp85_r1i1p1.ncml", season = c(12,1:11))
# EC.earth <- loadCMIP5(historical = "http://meteo.unican.es/tds5/dodsC/cmip5/EC-EARTH/EC-EARTH/historical/day/ec-earth_ec-earth_historical_r12i1p1.ncml",
#                       rcp8.5 = "http://meteo.unican.es/tds5/dodsC/cmip5/EC-EARTH/EC-EARTH/rcp85/day/ec-earth_ec-earth_rcp85_r12i1p1.ncml", season = c(12,1:11))
  ipsl <- loadCMIP5(historical = "/home/juan/Documents/2020_CMIP/Data/nc/psl_day_IPSL-CM5A-LR_historical_r1i1p1_19500101-20051231.nc",
                   rcp8.5 = "/home/juan/Documents/2020_CMIP/Data/nc/psl_day_IPSL-CM5A-LR_rcp85_r1i1p1_20060101-22051231.nc", season = c(12,1:11))
# mohc <- loadCMIP5(historical = "http://meteo.unican.es/tds5/dodsC/cmip5/MOHC/HADGEM2-ES/historical/day/mohc_hadgem2-es_historical_r1i1p1.ncml",
#                   rcp8.5 = "http://meteo.unican.es/tds5/dodsC/cmip5/MOHC/HADGEM2-ES/rcp85/day/mohc_hadgem2-es_rcp85_r1i1p1.ncml", season = c(12,1:11))
# mpi.esm.lr <- loadCMIP5(historical = "http://meteo.unican.es/tds5/dodsC/cmip5/MPI-M/MPI-ESM-LR/historical/day/mpi-m_mpi-esm-lr_historical_r1i1p1.ncml",
#                         rcp8.5 = "http://meteo.unican.es/tds5/dodsC/cmip5/MPI-M/MPI-ESM-LR/rcp85/day/mpi-m_mpi-esm-lr_rcp85_r1i1p1.ncml", season = c(12,1:11))
# mpi.esm.mr <- loadCMIP5(historical = "http://meteo.unican.es/tds5/dodsC/cmip5/MPI-M/MPI-ESM-MR/historical/day/mpi-m_mpi-esm-mr_historical_r1i1p1.ncml",
#                         rcp8.5 = "http://meteo.unican.es/tds5/dodsC/cmip5/MPI-M/MPI-ESM-MR/rcp85/day/mpi-m_mpi-esm-mr_rcp85_r1i1p1.ncml", season = c(12,1:11))
# ncc.nor <- loadCMIP5(historical = "http://meteo.unican.es/tds5/dodsC/cmip5/NCC/NORESM1-M/historical/day/ncc_noresm1-m_historical_r1i1p1.ncml",
#                      rcp8.5 = "http://meteo.unican.es/tds5/dodsC/cmip5/NCC/NORESM1-M/rcp85/day/ncc_noresm1-m_rcp85_r1i1p1.ncml", season = c(12,1:11))
# noaa.gfdl <- loadCMIP5(historical = "http://meteo.unican.es/tds5/dodsC/cmip5/NOAA-GFDL/GFDL-ESM2M/historical/day/noaa-gfdl_gfdl-esm2m_historical_r1i1p1.ncml",
#                        rcp8.5 = "http://meteo.unican.es/tds5/dodsC/cmip5/NOAA-GFDL/GFDL-ESM2M/rcp85/day/noaa-gfdl_gfdl-esm2m_rcp85_r1i1p1.ncml", season = c(12,1:11))
# 
  miroc <- loadGridData(dataset = "/home/juan/Documents/2020_CMIP/Data/CMIP5_MIROC_MIROC5_historical_r1i1p1.ncml",
              var = "psl",
              lonLim = lonLim,
              latLim = latLim,
              years = 1981:2010,
              season = c(12,1:11),
              time = "DD",
              aggr.d = "mean")
#
  grid1 <- loadGridData(dataset = "http://meteo.unican.es/tds5/dodsC/cmip5/MOHC/HADGEM2-ES/historical/day/mohc_hadgem2-es_historical_r1i1p1.ncml",
                        var = var,
                        lonLim = lonLim,
                        latLim = latLim,
                        years = 1980:2005,
                        time = "DD", 
                        aggr.d = "mean")
  
  grid2 <- loadGridData(dataset = "http://meteo.unican.es/tds5/dodsC/cmip5/MOHC/HADGEM2-ES/rcp85/day/mohc_hadgem2-es_rcp85_r1i1p1.ncml",
                        var = var,
                        lonLim = lonLim,
                        latLim = latLim,
                        years = 2005:2010,
                        time = "DD", 
                        aggr.d = "mean")
  hadgem <- bindGrid(grid1, grid2, dimension = "time")
  hadgem <- subsetGrid(hadgem, season = c(12,1:11))
 save(cccma, cnrm, EC.earth, ipsl, miroc, hadgem, mpi.esm.lr, ncc.nor, noaa.gfdl, file = "cmip5_europe.RData")
# 
# ### Fix grid to 10957 days. A different grid is needed to be loaded:
# grid <- cnrm
# if(getShape(cccma, dimension = "time") != 10957){
#   ind <- which((grid$Dates$start %in% cccma$Dates$start) == TRUE)
#   grid$Data[ind, , ] <- cccma$Data
#   cccma <- grid
# }
# 
# ### Lamb WTs from GCMs CMIP5: 
wts.cccma <- clusterGrid(cccma, type = "lamb")
wts.cnrm <- clusterGrid(cnrm, type = "lamb")
wts.EC.earth <- clusterGrid(EC.earth, type = "lamb")
wts.ipsl <- clusterGrid(ipsl, type = "lamb")
wts.miroc <- clusterGrid(miroc, type = "lamb")
wts.mohc <- clusterGrid(hadgem, type = "lamb")
wts.mpi.esm.lr <- clusterGrid(mpi.esm.lr, type = "lamb")
wts.ncc.nor <- clusterGrid(ncc.nor, type = "lamb")
wts.noaa.gfdl <- clusterGrid(noaa.gfdl, type = "lamb")

save(wts.cccma, wts.cnrm, wts.EC.earth, wts.ipsl, wts.miroc, wts.mohc, wts.mpi.esm.lr, wts.ncc.nor, wts.noaa.gfdl, file = "WTs_cmip5_europe.RData")

load("/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/WTs_cmip5.RData", verbose = TRUE)
load("/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/cmip5_europe.RData", verbose = TRUE)

##### Figure 3: Levelplot with RMSE between all the models (CMIP5, CMIP6, reanalysis products) and ERA-Interim

#### CMIP5 Relative Bias 

# getYearsAsINDEX(wts.cccma) %>% unique()
# a <- subsetGrid(wts.cccma,  season = c(12,1,2), years = getYearsAsINDEX(wts.cccma) %>% unique())
# getYearsAsINDEX(a) %>% range()
# getSeason(a)
# getSeason(wts.cccma)

# # #ERA-Interim Seasonal LWTs
# a <- freqWT(clusters.era.interim, season = DJF)
# b <- freqWT(clusters.era.interim, season = MAM)
# c <- freqWT(clusters.era.interim, season = JJA)
# d <- freqWT(clusters.era.interim, season = SON)
# 
#JRA Seasonal LWTs
a <- freqWT(clusters.jra, season = DJF)
b <- freqWT(clusters.jra, season = MAM)
c <- freqWT(clusters.jra, season = JJA)
d <- freqWT(clusters.jra, season = SON)

# #NCEP Seasonal LWTs
# a <- freqWT(clusters.ncep, season = DJF)
# b <- freqWT(clusters.ncep, season = MAM)
# c <- freqWT(clusters.ncep, season = JJA)
# d <- freqWT(clusters.ncep, season = SON)
# 
# #ERA-20C Seasonal LWTs
# a <- freqWT(clusters.era20, season = DJF)
# b <- freqWT(clusters.era20, season = MAM)
# c <- freqWT(clusters.era20, season = JJA)
# d <- freqWT(clusters.era20, season = SON)
      
### Frequencies Differences in LWTs

CMIP5.names <- c("CanESM2", "CNRM-CM5", "EC-EARTH", "GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5-LR", "MIROC5",
                 "MPI-ESM-LR", "NorESM1-M")
reanalysis <- c("ERA-Interim", "ERA-20C", "NCEP")

DJF.diff.freqs.cmip5 <- relativeFreqWT(clusters.era.interim, clusters.era20, clusters.ncep, wts.cccma, wts.cnrm, wts.EC.earth, wts.noaa.gfdl,
                                       wts.mohc, wts.ipsl, wts.miroc, wts.mpi.esm.lr, wts.ncc.nor, 
                                       ref = a,  cluster = WTs.index, season = DJF)

rownames(DJF.diff.freqs.cmip5) <- subsetWT
colnames(DJF.diff.freqs.cmip5) <- c(reanalysis, CMIP5.names)

MAM.diff.freqs.cmip5 <- relativeFreqWT(clusters.era.interim, clusters.era20, clusters.ncep, wts.cccma, wts.cnrm, wts.EC.earth, wts.noaa.gfdl,
                                       wts.mohc, wts.ipsl, wts.miroc, wts.mpi.esm.lr, wts.ncc.nor,
                                       ref = b, cluster = WTs.index, season = MAM)

rownames(MAM.diff.freqs.cmip5) <- subsetWT
colnames(DJF.diff.freqs.cmip5) <- c(reanalysis, CMIP5.names)

JJA.diff.freqs.cmip5 <- relativeFreqWT(clusters.era.interim, clusters.era20, clusters.ncep, wts.cccma, wts.cnrm, wts.EC.earth, wts.noaa.gfdl,
                                       wts.mohc, wts.ipsl, wts.miroc, wts.mpi.esm.lr, wts.ncc.nor,
                                       ref = c, cluster = WTs.index, season = JJA)

rownames(JJA.diff.freqs.cmip5) <- subsetWT
colnames(DJF.diff.freqs.cmip5) <- c(reanalysis, CMIP5.names)

SON.diff.freqs.cmip5 <- relativeFreqWT(clusters.era.interim, clusters.era20,clusters.ncep, wts.cccma, wts.cnrm, wts.EC.earth, wts.noaa.gfdl,
                                       wts.mohc, wts.ipsl, wts.miroc, wts.mpi.esm.lr, wts.ncc.nor,
                                       ref = d,cluster = WTs.index, season = SON)

rownames(SON.diff.freqs.cmip5) <- subsetWT
colnames(DJF.diff.freqs.cmip5) <- c(reanalysis, CMIP5.names)

yearly.diff.freqs.cmip5 <- relativeFreqWT(clusters.era.interim, clusters.era20,clusters.ncep, wts.cccma, wts.cnrm, wts.EC.earth, wts.noaa.gfdl,
                                       wts.mohc, wts.ipsl, wts.miroc, wts.mpi.esm.lr, wts.ncc.nor, 
                                       ref = freqWT(clusters.jra),  cluster = WTs.index)

rownames(yearly.diff.freqs.cmip5) <- subsetWT
colnames(DJF.diff.freqs.cmip5) <- c(reanalysis, CMIP5.names)

library(RColorBrewer)
library(gridExtra)

# display.brewer.all()
# RColorBrewer::brewer.pal(n = 9, "RdBu") %>% rev()
# pcolors <- RColorBrewer::brewer.pal(n = 9, "RdBu") %>% colorRampPalette()
# #  %>% print(split=c(2, 1, 2, 2), newpage=FALSE)
# p1 <- levelplot(DJF.diff.freqs.cmip5,
#           ylab = "CMIP5 GCM", xlab = "WTs", 
#           main = "(winter)",
#           col.regions = rev(pcolors(201)), 
#           at = seq(-15, 15, 0.5)) 
#  
# p2 <- levelplot(MAM.diff.freqs.cmip5,
#           ylab = "CMIP5 GCM", xlab = "WTs", 
#           main = "(spring)",
#           col.regions = rev(pcolors(201)), 
#           at = seq(-15, 15, 0.5)) 
# 
# p3 <- levelplot(JJA.diff.freqs.cmip5,
#           ylab = "CMIP5 GCM", xlab = "WTs", 
#           main = "(summer)",
#           col.regions = rev(pcolors(201)), 
#           at = seq(-15, 15, 0.5))  
# 
# p4 <- levelplot(SON.diff.freqs.cmip5,
#           ylab = "CMIP5 GCM", xlab = "WTs", 
#           main = "(fall)",
#           col.regions = rev(pcolors(201)), 
#           at = seq(-15, 15, 0.5))  
# #title(main = "Relative Freq. GCM vs ERA-Interim", outer = TRUE)
# 
# grid.arrange(p1, p2, p3, p4, ncol=2, top = "Relative Freq. GCM vs ERA-Interim")

bias.DJF.cmip5 <- matrix(nrow = nrow(DJF.diff.freqs.cmip5), 
                         ncol = ncol(DJF.diff.freqs.cmip5), 
                         dimnames = list(subsetWT, c(reanalysis, CMIP5.names)))
for(i in 1:ncol(DJF.diff.freqs.cmip5)){
  bias.DJF.cmip5[ ,i] <- DJF.diff.freqs.cmip5[,i]/a[subsetWT]
}

bias.MAM.cmip5 <- matrix(nrow = nrow(MAM.diff.freqs.cmip5), 
                         ncol = ncol(MAM.diff.freqs.cmip5), 
                         dimnames = list(subsetWT, c(reanalysis, CMIP5.names)))
for(i in 1:ncol(MAM.diff.freqs.cmip5)){
  bias.MAM.cmip5[ ,i] <- MAM.diff.freqs.cmip5[,i]/b[subsetWT]
}

bias.JJA.cmip5 <- matrix(nrow = nrow(JJA.diff.freqs.cmip5), 
                         ncol = ncol(JJA.diff.freqs.cmip5),
                         dimnames = list(subsetWT, c(reanalysis, CMIP5.names)))
for(i in 1:ncol(JJA.diff.freqs.cmip5)){
  bias.JJA.cmip5[ ,i] <- JJA.diff.freqs.cmip5[,i]/c[subsetWT]
}

bias.SON.cmip5 <- matrix(nrow = nrow(SON.diff.freqs.cmip5), 
                         ncol = ncol(SON.diff.freqs.cmip5), 
                         dimnames = list(subsetWT, c(reanalysis, CMIP5.names)))
for(i in 1:ncol(SON.diff.freqs.cmip5)){
  bias.SON.cmip5[ ,i] <- SON.diff.freqs.cmip5[,i]/d[subsetWT]
}

bias.yearly.cmip5 <- matrix(nrow = nrow(yearly.diff.freqs.cmip5), 
                         ncol = ncol(yearly.diff.freqs.cmip5), 
                         dimnames = list(subsetWT, c(reanalysis, CMIP5.names)))
for(i in 1:ncol(yearly.diff.freqs.cmip5)){
  bias.yearly.cmip5[ ,i] <- yearly.diff.freqs.cmip5[,i]/d[subsetWT]
}

#### CMIP6 Relative Bias  

# inst.model.run <- c("CCCma/CanESM5/historical/r1i1p1f1",
#                     "CNRM-CERFACS/CNRM-CM6-1/historical/r1i1p1f2",
#                     "IPSL/IPSL-CM6A-LR/historical/r1i1p1f1",
#                     "MIROC/MIROC6/historical/r1i1p1f1",
#                     "MOHC/HadGEM3-GC31-LL/historical/r1i1p1f3",
#                     "MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1",
#                     "NCC/NorESM2-LM/historical/r1i1p1f1")
# 
# dataset.cnrm <- "/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/CMIP6/CMIP6_CNRM-CERFACS_CNRM-CM6-1_historical_r1i1p1f2.ncml"
# dataset.cccma <- "/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/CMIP6/CMIP6_CCCma_CanESM5_historical_r1i1p1f1.ncml"
# dataset.ipsl <- "/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/CMIP6/CMIP6_IPSL_IPSL-CM6A-LR_historical_r1i1p1f1.ncml"
 dataset.miroc <- "https://data.meteo.unican.es/thredds/dodsC/devel/atlas/cmip6/CMIP/MIROC/MIROC-ES2L/historical/day/cmip6_CMIP_MIROC_MIROC-ES2L_historical_r1i1p1f2_day"
 "https://data.meteo.unican.es/thredds/dodsC/devel/atlas/cmip6/CMIP/MIROC/MIROC6/historical/day/cmip6_CMIP_MIROC_MIROC6_historical_r1i1p1f1_day"
# dataset.mohc <- "/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/CMIP6/CMIP6_MOHC_HadGEM3-GC31-LL_historical_r1i1p1f3.ncml"
# dataset.mpi <- "/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/CMIP6/CMIP6_MPI-M_MPI-ESM1-2-LR_historical_r1i1p1f1.ncml"
# dataset.ncc <- "/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/CMIP6/CMIP6_NCC_NorESM2-LM_historical_r1i1p1f1.ncml"
# dataset.gfdl <- "https://data.meteo.unican.es/thredds/dodsC/devel/atlas/cmip6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/day/cmip6_CMIP_NOAA-GFDL_GFDL-ESM4_historical_r1i1p1f1_day"
# dataset.ukesm1 <- "https://data.meteo.unican.es/thredds/dodsC/devel/atlas/cmip6/CMIP/MOHC/UKESM1-0-LL/historical/day/cmip6_CMIP_MOHC_UKESM1-0-LL_historical_r1i1p1f2_day"
# dataset.EC_earth3 <- "https://data.meteo.unican.es/thredds/dodsC/devel/atlas/cmip6/CMIP/EC-Earth-Consortium/EC-Earth3/historical/day/cmip6_CMIP_EC-Earth-Consortium_EC-Earth3_historical_r1i1p1f1_day"
#
# di <- dataInventory(dataset)
# 
# ipsl.cmip6 <- loadGridData(dataset = dataset.ipsl,
#                             var = "psl",
#                             lonLim = lonLim,
#                             latLim = latLim,
#                             years = 1981:2010,
#                             season = c(12,1:11),
#                             time = "DD", 
#                             aggr.d = "mean")
# 
miroc.cmip6 <- loadGridData(dataset = dataset.miroc,
                           var = "psl",
                           lonLim = lonLim,
                           latLim = latLim,
                           years = 1981:2010,
                           season = c(12,1:11),
                           time = "DD",
                           aggr.d = "mean")
# mohc.cmip6 <- loadGridData(dataset = dataset.mohc,
#                            var = "psl",
#                            lonLim = lonLim,
#                            latLim = latLim,
#                            years = 1981:2010,
#                            season = c(12,1:11),
#                            time = "DD", 
#                            aggr.d = "mean")
# mpi.cmip6 <- loadGridData(dataset = dataset.mpi,
#                            var = "psl",
#                            lonLim = lonLim,
#                            latLim = latLim,
#                            years = 1981:2010,
#                            season = c(12,1:11),
#                            time = "DD", 
#                            aggr.d = "mean")
# ncc.cmip6 <- loadGridData(dataset = dataset.ncc,
#                            var = "psl",
#                            lonLim = lonLim,
#                            latLim = latLim,
#                            years = 1981:2010,
#                            season = c(12,1:11),
#                            time = "DD", 
#                            aggr.d = "mean")
gfdl.cmip6 <- loadGridData(dataset = dataset.gfdl,
                           var = "psl",
                           lonLim = lonLim,
                           latLim = latLim,
                           years = 1981:2010,
                           season = c(12,1:11),
                           time = "DD", 
                           aggr.d = "mean")
ukesm1.cmip6 <- loadGridData(dataset = dataset.ukesm1,
                           var = "psl",
                           lonLim = lonLim,
                           latLim = latLim,
                           years = 1981:2010,
                           season = c(12,1:11),
                           time = "DD", 
                           aggr.d = "mean")
ec_earth.cmip6 <- loadGridData(dataset = dataset.EC_earth3,
                           var = "psl",
                           lonLim = lonLim,
                           latLim = latLim,
                           years = 1981:2010,
                           season = c(12,1:11),
                           time = "DD", 
                           aggr.d = "mean")
# 
# save(cccma.cmip6, cnrm.cmip6, ipsl.cmip6, miroc.cmip6, mohc.cmip6, mpi.cmip6, ncc.cmip6, file = "cmip6_europe.RData")
# 
# ### Clustering analysis: 
wts.cccma.cmip6 <- clusterGrid(cccma.cmip6, type = "lamb")
wts.cnrm.cmip6 <- clusterGrid(cnrm.cmip6, type = "lamb")
wts.ipsl.cmip6 <- clusterGrid(ipsl.cmip6, type = "lamb")
wts.miroc.cmip6 <- clusterGrid(miroc.cmip6, type = "lamb")
wts.ukesm1.cmip6 <- clusterGrid(ukesm1.cmip6, type = "lamb")
wts.mpi.lr.cmip6 <- clusterGrid(mpi.cmip6, type = "lamb")
wts.ncc.nor.cmip6 <- clusterGrid(ncc.cmip6, type = "lamb")
wts.ec_earth.cmip6.cmip6 <- clusterGrid(ec_earth.cmip6, type = "lamb")
wts.gfdl.cmip6 <- clusterGrid(gfdl.cmip6, type = "lamb")
# 
save(wts.cccma.cmip6, wts.cnrm.cmip6, wts.ec_earth.cmip6.cmip6, wts.gfdl.cmip6, wts.ipsl.cmip6, wts.miroc.cmip6, wts.ukesm1.cmip6, wts.mpi.lr.cmip6, wts.ncc.nor.cmip6, file = "wts_cmip6.RData")

load("/oceano/gmeteo/WORK/fdezja/2020_CMIP/Data/WTs_cmip6.RData", verbose = TRUE)


### Frequencies Differences in LWTs

WTs.index <- c(1,18,15,14,16,13,7,17)
subsetWT <- c("A", "C", "W", "SW", "NW", "S", "AW", "N")
CMIP6.names <- CMIP5.names <- c("CanESM2", "CNRM-CM5", "EC-EARTH", "GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5-LR", "MIROC5",
                                "MPI-ESM-LR", "NorESM1-M")

DJF.diff.freqs <- relativeFreqWT(wts.cccma.cmip6, wts.cnrm.cmip6, wts.ec_earth.cmip6.cmip6, wts.gfdl.cmip6, wts.ukesm1.cmip6, wts.ipsl.cmip6, wts.miroc.cmip6, 
                                 wts.mpi.lr.cmip6, wts.ncc.nor.cmip6, cluster = WTs.index, ref = a, season = DJF)

rownames(DJF.diff.freqs) <- subsetWT
colnames(DJF.diff.freqs) <- CMIP6.names

MAM.diff.freqs <- relativeFreqWT(wts.cccma.cmip6, wts.cnrm.cmip6, wts.ec_earth.cmip6.cmip6, wts.gfdl.cmip6, wts.ukesm1.cmip6, wts.ipsl.cmip6, wts.miroc.cmip6, 
                                 wts.mpi.lr.cmip6, wts.ncc.nor.cmip6, cluster = WTs.index, ref = b, season = MAM)

rownames(MAM.diff.freqs) <- subsetWT
colnames(MAM.diff.freqs) <- CMIP6.names

JJA.diff.freqs <- relativeFreqWT(wts.cccma.cmip6, wts.cnrm.cmip6, wts.ec_earth.cmip6.cmip6, wts.gfdl.cmip6, wts.ukesm1.cmip6, wts.ipsl.cmip6, wts.miroc.cmip6, 
                                 wts.mpi.lr.cmip6, wts.ncc.nor.cmip6, cluster = WTs.index, ref = c, season = JJA)

rownames(JJA.diff.freqs) <- subsetWT
colnames(JJA.diff.freqs) <- CMIP6.names

SON.diff.freqs <- relativeFreqWT(wts.cccma.cmip6, wts.cnrm.cmip6, wts.ec_earth.cmip6.cmip6, wts.gfdl.cmip6, wts.ukesm1.cmip6, wts.ipsl.cmip6, wts.miroc.cmip6, 
                                 wts.mpi.lr.cmip6, wts.ncc.nor.cmip6, ref = d, cluster = WTs.index, season = SON)

rownames(SON.diff.freqs) <- subsetWT
colnames(SON.diff.freqs) <- CMIP6.names

yearly.diff.freqs <- relativeFreqWT(wts.cccma.cmip6, wts.cnrm.cmip6, wts.ec_earth.cmip6.cmip6, wts.gfdl.cmip6, wts.ukesm1.cmip6, wts.ipsl.cmip6, wts.miroc.cmip6, 
                                    wts.mpi.lr.cmip6, wts.ncc.nor.cmip6, ref = freqWT(clusters.era.interim),  cluster = WTs.index)

rownames(yearly.diff.freqs) <- subsetWT
colnames(yearly.diff.freqs) <- CMIP6.names

# RColorBrewer::brewer.pal(n = 9, "RdBu") %>% rev()
# pcolors <- RColorBrewer::brewer.pal(n = 9, "RdBu") %>% colorRampPalette()
# p1 <- levelplot(DJF.diff.freqs,
#                 ylab = "CMIP6 GCM", xlab = "WTs", 
#                 main = "(winter)",
#                 col.regions = rev(pcolors(201)), 
#                 at = seq(-15, 15, 0.5)) 
# 
# p2 <- levelplot(MAM.diff.freqs,
#                 ylab = "CMIP6 GCM", xlab = "WTs", 
#                 main = "(spring)",
#                 col.regions = rev(pcolors(201)), 
#                 at = seq(-15, 15, 0.5))  
# 
# p3 <- levelplot(JJA.diff.freqs,
#                 ylab = "CMIP6 GCM", xlab = "WTs", 
#                 main = "(summer)",
#                 col.regions = rev(pcolors(201)), 
#                 at = seq(-15, 15, 0.5))   
# 
# p4 <- levelplot(SON.diff.freqs,
#                 ylab = "CMIP6 GCM", xlab = "WTs", 
#                 main = "(fall)",
#                 col.regions = rev(pcolors(201)), 
#                 at = seq(-15, 15, 0.5))   
# #title(main = "Relative Freq. GCM vs ERA-Interim", outer = TRUE)
# 
# grid.arrange(p1, p2, p3, p4, ncol=2, top = "Relative Freq. GCM vs ERA-Interim")

### Relative Bias CMIP6:

bias.DJF.cmip6 <- matrix(nrow = nrow(DJF.diff.freqs), ncol = ncol(DJF.diff.freqs), dimnames = list(subsetWT, CMIP6.names))
for(i in 1:ncol(DJF.diff.freqs)){
  bias.DJF.cmip6[ ,i] <- DJF.diff.freqs[,i]/a[subsetWT]
}

bias.MAM.cmip6 <- matrix(nrow = nrow(MAM.diff.freqs), ncol = ncol(MAM.diff.freqs), dimnames = list(subsetWT, CMIP6.names))
for(i in 1:ncol(MAM.diff.freqs)){
  bias.MAM.cmip6[ ,i] <- MAM.diff.freqs[,i]/b[subsetWT]
}

bias.JJA.cmip6 <- matrix(nrow = nrow(JJA.diff.freqs), ncol = ncol(JJA.diff.freqs), dimnames = list(subsetWT, CMIP6.names))
for(i in 1:ncol(JJA.diff.freqs)){
  bias.JJA.cmip6[ ,i] <- JJA.diff.freqs[,i]/c[subsetWT]
}

bias.SON.cmip6 <- matrix(nrow = nrow(SON.diff.freqs), ncol = ncol(SON.diff.freqs), dimnames = list(subsetWT, CMIP6.names))
for(i in 1:ncol(SON.diff.freqs)){
  bias.SON.cmip6[ ,i] <- SON.diff.freqs[,i]/d[subsetWT]
}

bias.yearly.cmip6 <- matrix(nrow = nrow(yearly.diff.freqs), 
                            ncol = ncol(yearly.diff.freqs), 
                            dimnames = list(subsetWT, CMIP6.names))
for(i in 1:ncol(yearly.diff.freqs)){
  bias.yearly.cmip6[ ,i] <- yearly.diff.freqs[,i]/d[subsetWT]
}

## Relative Bias Data-Frame: 

order.season <- c("DJF","MAM","JJA","SON", "Year")
GCM_evaluation.DF <- data.frame(WT= factor(rep(subsetWT, 21*5),levels = subsetWT),
                                GCM = factor(c(rep(c(rep("NorESM1-M", 8),
                                               rep("MPI-ESM-LR", 8), 
                                               rep("MIROC5", 8), 
                                               rep("IPSL-CM5-LR", 8),
                                               rep("HadGEM2-ES", 8), 
                                               rep("GFDL-ESM2M", 8),
                                               rep("EC-EARTH", 8), 
                                               rep("CNRM-CM5", 8), 
                                               rep("CanESM2", 8)),5*2),
                                             rep(c(rep("NCEP", 8),
                                               rep("ERA-20C", 8),
                                               rep("JRA", 8)),5)), levels = c(reanalysis, CMIP5.names[9:1])), 
                                Experiment = factor(c(rep("CMIP6", 9*8*5), rep("CMIP5", 12*8*5)),levels = c("CMIP6","CMIP5")),
                                Season = factor(c(rep("DJF",8*9), rep("MAM",8*9), rep("JJA",8*9), rep("SON",8*9), rep("Yearly",8*9),
                                           rep("DJF",8*12), rep("MAM",8*12), rep("JJA",8*12), rep("SON",8*12), rep("Yearly",8*12)),levels = order.season),
                                Bias_rel = c(bias.DJF.cmip6[ ,9:1], bias.MAM.cmip6[ ,9:1], 
                                             bias.JJA.cmip6[,9:1], bias.SON.cmip6[,9:1], bias.yearly.cmip6[ ,9:1],
                                             bias.DJF.cmip5[,12:1], bias.MAM.cmip5[,12:1], 
                                             bias.JJA.cmip5[,12:1], bias.SON.cmip5[,12:1],bias.yearly.cmip5[ ,12:1])) 

write.csv(GCM_evaluation.DF, file = "GCM_evaluation.csv")
save(GCM_evaluation.DF, file = "GCM_evaluation.DF.RData")

# pcolors <- RColorBrewer::brewer.pal(n = 9, "RdBu") %>% colorRampPalette()
# 
# lattice::levelplot(Bias_rel ~ WT * GCM | Season * Experiment , data = GCM_evaluation.DF,
#                    main = "Relative Bias - GCM vs Reanalysis",
#                    col.regions = pcolors(201) %>% rev(),
#                    at = seq(-1.4, 1.4, 0.1), 
#                    scales=list(alternating=1),
#                    par.settings=list(panel.background = list(col = "grey85")))


### The mean-square Error (RMSE): 

GCM.cmip5.rmse <- matrix(ncol = 5, nrow = 12, dimnames = list(c(reanalysis, CMIP5.names), order.season))
GCM.cmip6.rmse <- matrix(ncol = 5, nrow = 9, dimnames = list(CMIP6.names, order.season))

for (i in 1:nrow(GCM.cmip5.rmse)) {
  GCM.cmip5.rmse[i,1] <- rmse(bias.DJF.cmip5[,i])
  GCM.cmip5.rmse[i,2] <- rmse(bias.MAM.cmip5[ ,i])
  GCM.cmip5.rmse[i,3] <- rmse(bias.JJA.cmip5[ ,i])
  GCM.cmip5.rmse[i,4] <- rmse(bias.SON.cmip5[ ,i])
  GCM.cmip5.rmse[i,5] <- rmse(bias.yearly.cmip5[ ,i])
}
for (i in 1:nrow(GCM.cmip6.rmse)) {
  GCM.cmip6.rmse[i,1] <- rmse(bias.DJF.cmip6[ ,i])
  GCM.cmip6.rmse[i,2] <- rmse(bias.MAM.cmip6[ ,i])
  GCM.cmip6.rmse[i,3] <- rmse(bias.JJA.cmip6[ ,i])
  GCM.cmip6.rmse[i,4] <- rmse(bias.SON.cmip6[ ,i])
  GCM.cmip6.rmse[i,5] <- rmse(bias.yearly.cmip6[ ,i])
}

#Create GCM.rmse.DF:
GCM.rmse.DF <- data.frame(GCM = factor(c(rep(reanalysis, 5), rep(CMIP5.names,5) , rep(CMIP6.names,5)), levels = c(reanalysis, CMIP5.names[10:1])),
                         Experiment = factor(c(rep("Reanalysis", 3*5), rep("CMIP5", 9*5), rep("CMIP6", 9*5)),levels = c("Reanalysis", "CMIP5","CMIP6")),
                         Season = factor(c(rep("DJF",3),rep("MAM",3),rep("JJA",3),rep("SON",3), rep("Year",3),
                                           rep("DJF",9),rep("MAM",9),rep("JJA",9),rep("SON",9), rep("Year",9),
                                           rep("DJF",9),rep("MAM",9),rep("JJA",9),rep("SON",9), rep("Year",9)), levels = order.season),
                         RMSE = c(as.vector(GCM.cmip5.rmse[1:3, ]), as.vector(GCM.cmip5.rmse[4:12, ]), as.vector(GCM.cmip6.rmse)))

save(GCM.rmse.DF, file = "GCM.rmse.RData")
write.csv(GCM.rmse.DF, file = "GCM.rmse.csv")

dev.new()
pcolors <- RColorBrewer::brewer.pal(n = 9, "OrRd") %>% colorRampPalette()


z <- subset(GCM.rmse.DF, Experiment == "Reanalysis")
p1 <- lattice::levelplot(RMSE ~ Season + GCM , data = z,
                         col.regions = pcolors(201),
                         at = seq(0, 0.6, 0.02),
                         colorkey = list(space = 'bottom'),
                         par.settings=list(layout.widths=list(key.ylab.padding = 2))) +  
      lattice::xyplot(RMSE ~ Season + GCM , data = z,
                      panel = function(x, y, ...) {
                      ltext(x = z$Season, y = z$GCM, labels = round(z$RMSE, 2), cex = 0.9, font = 1,
                      fontfamily = "HersheySans")})

w <- subset(GCM.rmse.DF, Experiment == "CMIP5")
p2 <- lattice::levelplot(RMSE ~ Season + GCM , data = w,
                         col.regions = pcolors(201),
                         at = seq(0, 0.6, 0.02),
                         scales=list(alternating=1),
                         colorkey = list(space = 'bottom')) +  
  lattice::xyplot(RMSE ~ Season + GCM , data = w,
                  panel = function(x, y, ...) {
                    ltext(x = w$Season, y = factor(CMIP5.names, levels = CMIP5.names[9:1]), labels = round(w$RMSE, 2), cex = 0.9, font = 1,
                          fontfamily = "HersheySans")})

v <- subset(GCM.rmse.DF, Experiment == "CMIP6")
p3 <- lattice::levelplot(RMSE ~ Season + GCM , data = v,
                         col.regions = pcolors(201),
                         at = seq(0, 0.6, 0.02),
                         scales=list(y = list(alternating=2)),
                         colorkey = list(space = 'bottom')) +  
  lattice::xyplot(RMSE ~ Season + GCM , data = v,
                  panel = function(x, y, ...) {
                    ltext(x = v$Season, y = factor(CMIP5.names, levels = CMIP5.names[9:1]), labels = round(v$RMSE, 2), cex = 0.9, font = 1,
                          fontfamily = "HersheySans")})

comb_levObj <- c(Reanalysis = p1, CMIP5 = p2, CMIP6 = p3, layout = c(3, 1))
print(comb_levObj)
update(comb_levObj, main = "Root Mean Square Error : GCM & Reanalysis vs ERA-Interim", xlab = NULL,
       scales = list(x = list(alternating=3), y = list(rot=0)))

  
  # lattice::levelplot(RMSE ~ Season + GCM | Experiment, data = GCM.rmse.DF, 
#                     main = "Root Mean Square Error - GCM vs Reanalysis",
#                     col.regions = pcolo rs(201),
#                     at = seq(0, 0.6, 0.02),
#                     scales=list(alternating=1),
#                     layout = c(3,1), 
#                     labels = FALSE,
#                     par.settings=list(panel.background = list(col = "grey85")))


###### FIGURE 4: Relative Bias 8 ordered WTs: -------------------------------------
#Order GCMs:

indices <- matrix(nrow = 21, ncol = 5, dimnames = list(NULL,order.season))
# year.order <- vector("numeric", length = 21)
# DJF.order <- vector("numeric", length = 21)
# MAM.order <- vector("numeric", length = 21)
# JJA.order <- vector("numeric", length = 21)
# SON.order <- vector("numeric", length = 21)
KL.yearly <- subset(KLdivergence.DF, Season == "Year")
KL.yearly.order <- KL.yearly[order(KL.yearly$KL),]
indices[,5] <- as.numeric(row.names(KL.yearly.order))
KL.DJF <- subset(KLdivergence.DF, Season == "DJF")
KL.DJF.order <- KL.DJF[order(KL.DJF$KL),]
indices[,1] <- as.numeric(row.names(KL.DJF.order))
KL.MAM <- subset(KLdivergence.DF, Season == "MAM")
KL.MAM.order <- KL.MAM[order(KL.MAM$KL),]
indices[,2] <- as.numeric(row.names(KL.MAM.order))
KL.JJA <- subset(KLdivergence.DF, Season == "JJA")
KL.JJA.order <- KL.JJA[order(KL.JJA$KL),]
indices[,3] <- as.numeric(row.names(KL.JJA.order))
KL.SON <- subset(KLdivergence.DF, Season == "SON")
KL.SON.order <- KL.SON[order(KL.SON$KL),]
indices[,4] <- as.numeric(row.names(KL.SON.order))

indices.aux <- indices

for(i in 1:length(order.season)){
  for (j in 1:nrow(indices)) {
    if(i == 1){ #Es DJF
      if((indices[j,i] < 70 && indices[j,i] > 60)){ #es CMIP6
        indices.aux[j,i] <- as.numeric(indices[j,i]) - 48
      }
      if((indices[j,i] < 25 && indices[j,i] > 15)){ #es CMIP5
        indices.aux[j,i] <- indices[j,i] - 15 + 3
      }
      if((indices[j,i] < 4 && indices[j,i] > 0)){ #es reanálisis
        indices.aux[j,i] <- indices[j,i]
      }
    }
    if(i == 2){ #Es MAM
      if((indices[j,i] < 79 && indices[j,i] > 69)){ #es CMIP6
        indices.aux[j,i] <- indices[j,i] - 69 + 12
      }
      if((indices[j,i] < 34 && indices[j,i] > 24)){ #es CMIP5
        indices.aux[j,i] <- indices[j,i] - 24 + 3
      }
      if((indices[j,i] < 7 && indices[j,i] > 3)){ #es reanálisis
        indices.aux[j,i] <- indices[j,i] - 3
      }
    }
    if(i == 3){ #Es JJA
      if((indices[j,i] < 88 && indices[j,i] > 78)){ #es CMIP6
        indices.aux[j,i] <- indices[j,i] - 78 + 12
      }
      if((indices[j,i] < 43 && indices[j,i] > 33)){ #es CMIP5
        indices.aux[j,i] <- indices[j,i] - 33 + 3
      }
      if((indices[j,i] < 10 && indices[j,i] > 6)){ #es reanálisis
        indices.aux[j,i] <- indices[j,i] -6
      }
    }
    if(i == 4){ #Es SON
      if((indices[j,i] < 97 && indices[j,i] > 87)){ #es CMIP6
        indices.aux[j,i] <- indices[j,i] - 87 + 12
      }
      if((indices[j,i] < 52 && indices[j,i] > 42)){ #es CMIP5
        indices.aux[j,i] <- indices[j,i] - 42 + 3
      }
      if((indices[j,i] < 13 && indices[j,i] > 9)){ #es reanálisis
        indices.aux[j,i] <- indices[j,i] - 9
      }
    }
    if(i == 5){ #Es year
      if((indices[j,i] < 106 && indices[j,i] > 96)){ #es CMIP6
        indices.aux[j,i] <- indices[j,i] - 96 + 12
      }
      if((indices[j,i] < 61 && indices[j,i] > 51)){ #es CMIP5
        indices.aux[j,i] <- indices[j,i] - 51 + 3
      }
      if((indices[j,i] < 16 && indices[j,i] > 12)){ #es reanálisis
        indices.aux[j,i] <- indices[j,i] - 12
      }
    }
  }
}

KL.yearly <- subset(KLdivergence.DF, Season == "Year")
KL.yearly.order <- KL.yearly[order(KL.yearly$KL),]
index <- as.numeric(row.names(KL.yearly.order))

KL.DJF <- subset(KLdivergence.DF, Season == "DJF")
KL.DJF.order <- KL.DJF[order(KL.yearly$KL),]
index.DJF <- as.numeric(row.names(KL.DJF.order))
KL.MAM <- subset(KLdivergence.DF, Season == "MAM")
KL.MAM.order <- KL.MAM[order(KL.yearly$KL),]
index.MAM <- as.numeric(row.names(KL.MAM.order))
KL.JJA <- subset(KLdivergence.DF, Season == "JJA")
KL.JJA.order <- KL.JJA[order(KL.yearly$KL),]
index.JJA <- as.numeric(row.names(KL.JJA.order))
KL.SON <- subset(KLdivergence.DF, Season == "SON")
KL.SON.order <- KL.SON[order(KL.yearly$KL),]
index.SON <- as.numeric(row.names(KL.SON.order))

GCM_evaluation.mat <- cbind(bias.DJF.cmip5[,1:3], bias.MAM.cmip5[,1:3], bias.JJA.cmip5[,1:3], bias.SON.cmip5[,1:3], bias.yearly.cmip5[,1:3], 
                            bias.DJF.cmip5[,4:12], bias.MAM.cmip5[,4:12], bias.JJA.cmip5[,4:12], bias.SON.cmip5[,4:12], bias.yearly.cmip5[,4:12],
                            bias.DJF.cmip6[,], bias.MAM.cmip6[,], bias.JJA.cmip6[,], bias.SON.cmip6[,], bias.yearly.cmip6[,])
pcolors <- RColorBrewer::brewer.pal(n = 9, "RdBu") %>% colorRampPalette()

### Compute proportion Z-test among GCM-WTs' Freqs and Reanalysis WTs' Freq 

names.cmip5 <- c("CanESM2", "CNRM-CM5", "EC-EARTH", "GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5-LR", "MIROC5", "MPI-ESM-LR", "NorESM1-M")
names.cmip6 <- c("CanESM5", "CNRM-CM6-1","EC-EARTH3", "GFDL-ESM4", "UKESM1-0-LL", "IPSL-CM6A-LR", "MIROC6",  "MPI-ESM1-2-LR", "NorESM2-LM")
GCM.names.year.order <- c(reanalysis, names.cmip5, names.cmip6)[indices.aux[,5]]

#CROSS seasonal order with year order: 
names.cmip5.2 <- c("CanESM2", "CNRM-CM5", "EC-EARTH", "GFDL-ESM2M", "HadGEM2-ES", "IPSL-CM5-LR", "MIROC5", "MPI-ESM-LR", "NorESM1-M")
names.cmip6.2 <- c("CanESM5", "CNRM-CM6-1","EC-EARTH3", "GFDL-ESM4", "UKESM1-0-LL", "IPSL-CM6A-LR", "MIROC6",  "MPI-ESM1-2-LR", "NorESM2-LM")
GCM.names.2 <- c(reanalysis, names.cmip5.2, names.cmip6.2)

GCM.names <- GCM.names.year.order
aux <- match(GCM.names.year.order, GCM.names.2[indices.aux[,1]])
GCM.names <- paste0(GCM.names, " (", aux, ",")
aux <- match(GCM.names.year.order, GCM.names.2[indices.aux[,2]])
GCM.names <- paste0(GCM.names, " ", aux, ",")
aux <- match(GCM.names.year.order, GCM.names.2[indices.aux[,3]])
GCM.names <- paste0(GCM.names, " ", aux, ",")
aux <- match(GCM.names.year.order, GCM.names.2[indices.aux[,4]])
GCM.names <- paste0(GCM.names, " ", aux, ")")

data <- matrix(t(GCM_evaluation.mat)[index.DJF, ], 
               ncol = ncol(t(GCM_evaluation.mat)[index.DJF, ]), 
               nrow = nrow(t(GCM_evaluation.mat)[index.DJF, ]), 
               dimnames = list(GCM.names.year.order, subsetWT)) %>% t()
colors.y.p1 <- c(rep("black",2), "black", "blue","red","blue","red","blue","blue","blue","blue","blue","red","red","red","red","blue","red","red","blue","red")
points.p1 <- freqWT.pvalues(data.mat = data, season = DJF, clusters = subsetWT, ref = clusters.jra, alpha = 0.05)
points.p1[ ,2] <- abs(points.p1[ ,2]-22)
colnames(data) <- GCM.names
colorkey.labels <- c(-1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
p1 <- lattice::levelplot(data[ ,21:1], 
                         xlab = "WTs", ylab = "GCM",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1),
                         colorkey = list(space = 'bottom',
                                         labels=list(at = seq(-1.2, 1.2, 0.2), 
                                                     labels = colorkey.labels)),
                         par.settings = list(layout.heights = list(xlab.key.padding = 2)),
                         scales = list(y = list(col = rev(colors.y.p1))) ) +
  layer(sp::sp.points(sp::SpatialPoints(points.p1), pch= 4, col = 1, cex=2))

data <- matrix(t(GCM_evaluation.mat)[index.MAM, ], 
               ncol = ncol(t(GCM_evaluation.mat)[index.MAM, ]), 
               nrow = nrow(t(GCM_evaluation.mat)[index.MAM, ]), 
               dimnames = list(GCM.names.year.order, subsetWT)) %>% t()
points.p2 <- freqWT.pvalues(data.mat = data, season = MAM, clusters = subsetWT, ref = clusters.jra, alpha = 0.05)
points.p2[ ,2] <- abs(points.p2[ ,2]-22)
colnames(data) <- GCM.names

p2 <- lattice::levelplot(data[ ,21:1], 
                         xlab = "WTs", ylab = "GCM",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1), 
                         colorkey = list(space = 'left')) +
  layer(sp::sp.points(sp::SpatialPoints(points.p2), pch= 4, col = 1, cex=2))

data <- matrix(t(GCM_evaluation.mat)[index.JJA, ], 
               ncol = ncol(t(GCM_evaluation.mat)[index.JJA, ]), 
               nrow = nrow(t(GCM_evaluation.mat)[index.JJA, ]), 
               dimnames = list(GCM.names.year.order, subsetWT)) %>% t()
points.p3 <- freqWT.pvalues(data.mat = data, season = JJA, clusters = subsetWT, ref = clusters.jra, alpha = 0.05)
points.p3[ ,2] <- abs(points.p3[ ,2]-22)
colnames(data) <- GCM.names

p3 <- lattice::levelplot(data[ ,21:1], 
                         xlab = "WTs", ylab = "GCM",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1), 
                         colorkey = list(space = 'left')) +
  layer(sp::sp.points(sp::SpatialPoints(points.p3), pch= 4, col = 1, cex=2))

data <- matrix(t(GCM_evaluation.mat)[index.SON, ], 
               ncol = ncol(t(GCM_evaluation.mat)[index.SON, ]), 
               nrow = nrow(t(GCM_evaluation.mat)[index.SON, ]), 
               dimnames = list(GCM.names.year.order, subsetWT)) %>% t()
points.p4 <- freqWT.pvalues(data.mat = data, season = SON, clusters = subsetWT, ref = clusters.jra, alpha = 0.05)
points.p4[ ,2] <- abs(points.p4[ ,2]-22)
colnames(data) <- GCM.names

p4 <- lattice::levelplot(data[ ,21:1], 
                         xlab = "WTs", ylab = "GCM",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1),
                         colorkey = list(space = 'left')) +
  layer(sp::sp.points(sp::SpatialPoints(points.p4), pch= 4, col = 1, cex=2))

comb_levObj <- c(DJF = p1, MAM = p2, JJA = p3, SON = p4, layout = c(4, 1))
print(comb_levObj)
update(comb_levObj, main = "Rel. Bias GCM & Reanalysis vs JRA",
       scales = list(tck = c(1,0), 
                     x = list(alternating=1), 
                     y = list(rot=0, col = rev(colors.y.p1))))



#Seasonal ordered:
GCM.rmse.DF.ordered <- GCM.rmse.DF[order(GCM.rmse.DF$RMSE),]
print(GCM.rmse.DF.ordered)
index <- as.numeric(row.names(GCM.rmse.DF.ordered))

GCM_evaluation.mat <- cbind(bias.DJF.cmip5[,1:3], bias.MAM.cmip5[,1:3], bias.JJA.cmip5[,1:3], bias.SON.cmip5[,1:3], 
                            bias.DJF.cmip5[,4:13], bias.MAM.cmip5[,4:13], bias.JJA.cmip5[,4:13], bias.SON.cmip5[,4:13], 
                            bias.DJF.cmip6[,], bias.MAM.cmip6[,], bias.JJA.cmip6[,], bias.SON.cmip6[,])
# GCM_evaluation.mat <- cbind(DJF.diff.freqs.cmip5[,], MAM.diff.freqs.cmip5[,], JJA.diff.freqs.cmip5[,], SON.diff.freqs.cmip5[,], 
#                             DJF.diff.freqs[,], MAM.diff.freqs[,], JJA.diff.freqs[,], SON.diff.freqs[,])

dim(GCM_evaluation.mat)

GCM.ordered <- GCM_evaluation.mat[ ,index] %>%  t()
names.GCM.ordered <- paste0(row.names(t(GCM_evaluation.mat[ ,index])), " - ", GCM.rmse.DF.ordered$Season, " - ", GCM.rmse.DF.ordered$Experiment)

row.names(GCM.ordered) <- names.GCM.ordered
pcolors <- RColorBrewer::brewer.pal(n = 9, "RdBu") %>% colorRampPalette()
lattice::levelplot(t(GCM.ordered), 
                   xlab = "WTs", ylab = "GCM - Season - Experiment",
                   main = "Bias GCM vs Reanalysis",
                   col.regions = rev(pcolors(201)), 
                   at = seq(-1.4, 1.4, 0.1))
#at = seq(-3.4, 3.4, 0.1)

##Filter rows per season in last data matrix:
#DJF: 

GCM.DJF <- grep("DJF", names.GCM.ordered)
GCM.MAM <- grep("MAM", names.GCM.ordered)
GCM.JJA <- grep("JJA", names.GCM.ordered)
GCM.SON <- grep("SON", names.GCM.ordered)

p1 <- lattice::levelplot(t(GCM.ordered[GCM.DJF, ]), 
                         xlab = "WTs", ylab = "GCM - DJF - Experiment",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1))
p2 <- lattice::levelplot(t(GCM.ordered[GCM.MAM, ]), 
                         xlab = "WTs", ylab = "GCM - MAM - Experiment",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1))
p3 <- lattice::levelplot(t(GCM.ordered[GCM.JJA, ]), 
                         xlab = "WTs", ylab = "GCM - JJA - Experiment",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1))
p4 <- lattice::levelplot(t(GCM.ordered[GCM.SON, ]), 
                         scales=list(x=list(cex=0.8),y=list(cex=0.8)),
                         xlab = "WTs", ylab = "GCM - SON - Experiment",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1))

gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2, top = "Bias Rel. GCM vs Reanalysis")


# Compute proportion Z-test among GCM-WTs' Freqs and Reanalysis WTs' Freq 

names.cmip5.2 <- c("CanESM2", "CNRM-CM5", "EC-EARTH", "IPSL-CM5A-MR", "MIROC5",
                 "HadGEM2-ES", "MPI-ESM-LR", "MPI-ESM-MR", "NorESM1-M", "GFDL-ESM2M")
names.cmip6.2 <- c("CanESM5", "CNRM-CM6-1", "IPSL-CM6A-LR", "MIROC6","HadGEM3-GC31-LL", "MPI-ESM1-2-LR", "NorESM2-LM")
GCM.names.2 <- c("JRA", "ERA-20C", "NCEP", names.cmip5.2, names.cmip6.2)

aux <- subset(GCM.rmse.DF.ordered, Season == "DJF")
data <- matrix(t(GCM.ordered[GCM.DJF, ]), ncol = 20, nrow = 8, 
               dimnames = list( subsetWT, GCM.names.2[DJF.order]))
colors.y.p1 <- c(rep("black",2), "blue","blue","red", "black","red","blue","blue","blue","red","blue","blue","blue","blue","red","red","red","red","red","red")
points.p1 <- freqWT.pvalues(data.mat = data, season = DJF, clusters = subsetWT, ref = clusters.era.interim, alpha = 0.05)
points.p1[ ,2] <- abs(points.p1[ ,2]-22)

p1 <- lattice::levelplot(data[ ,20:1], 
                         xlab = "WTs", ylab = "GCM",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1),
                         colorkey = list(space = 'left'),
                         par.settings=list(layout.widths=list(key.ylab.padding = 2)),
                         scales=list(y=list(col = rev(colors.y.p1))) 
                         ) +
  layer(sp.points(SpatialPoints(points.p1), pch= 4, col = 1, cex=2))

aux <- subset(GCM.rmse.DF.ordered, Season == "MAM")
data <- matrix(t(GCM.ordered[GCM.MAM, ]), ncol = 20, nrow = 8, 
               dimnames = list( subsetWT, GCM.names.2[MAM.order]))
colors.y.p2 <- c("blue","blue","red","blue","red","blue","red","red","red","blue","blue","blue","red","red","red","red","red")
points.p2 <- freqWT.pvalues(data.mat = data, season = MAM, clusters = subsetWT, ref = clusters.era.interim, alpha = 0.05)
points.p2[ ,2] <- abs(points.p2[ ,2]-21)
p2 <- lattice::levelplot(data[ ,20:1], 
                         xlab = "WTs", ylab = "GCM",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1), 
                         colorkey = list(space = 'left'),
                         scales=list(y=list(col = colors.y.p2)), 
                         ) +
  layer(sp.points(SpatialPoints(points.p2), pch= 4, col = 1, cex=2))

aux <- subset(GCM.rmse.DF.ordered, Season == "JJA")
data <- matrix(t(GCM.ordered[GCM.JJA, ]), ncol = 20, nrow = 8, 
               dimnames = list( subsetWT, GCM.names.2[JJA.order]))
colors.y.p3 <- c("red","red","blue","red","red","red","blue","blue","blue","red","red","blue","red","red","red","blue","blue")
points.p3 <- freqWT.pvalues(data.mat = data, season = JJA, clusters = subsetWT, ref = clusters.era.interim, alpha = 0.05)
points.p3[ ,2] <- abs(points.p3[ ,2]-21)

p3 <- lattice::levelplot(data[ ,20:1], 
                         xlab = "WTs", ylab = "GCM",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1), 
                         colorkey = list(space = 'left'),
                         scales=list(y=list(col = rev(colors.y.p3))), 
                         ) +
  layer(sp.points(SpatialPoints(points.p3), pch= 4, col = 1, cex=2))

aux <- subset(GCM.rmse.DF.ordered, Season == "SON")
data <- matrix(t(GCM.ordered[GCM.SON, ]), ncol = 20, nrow = 8, 
               dimnames = list( subsetWT, GCM.names.2[SON.order]))
colors.y.p4 <- c("blue","blue","red","red","red","red","blue","blue","blue","red","red","red","red","blue","red","blue","red")
points.p4 <- freqWT.pvalues(data.mat = data, season = SON, clusters = subsetWT, ref = clusters.era.interim, alpha = 0.05)
points.p4[ ,2] <- abs(points.p4[ ,2]-21)
p4 <- lattice::levelplot(data[ ,20:1], 
                         xlab = "WTs", ylab = "GCM",
                         col.regions = rev(pcolors(201)), 
                         at = seq(-1.4, 1.4, 0.1),
                         colorkey = list(space = 'left'),
                         scales=list(y=list(col = rev(colors.y.p4))), 
                         ) +
  layer(sp.points(SpatialPoints(points.p4), pch= 4, col = 1, cex=2))

#print(p1+p2+p3+p4)
comb_levObj <- c(JJA = p3, SON = p4, DJF = p1, MAM = p2, layout = c(2, 2))
print(comb_levObj)
update(comb_levObj, main = "Rel. Bias GCM vs Reanalysis",
       scales = list(x = list(alternating=3), y = list(rot=0, col = c(rev(colors.y.p1), colors.y.p2, rev(colors.y.p3), rev(colors.y.p4)))))



  ## Matriz de transición -----------------------------------------
source('/oceano/gmeteo/WORK/fdezja/2020_CMIP/Scripts/transitionProb.R')


# Plotting 
pcolors <- RColorBrewer::brewer.pal(n = 9, "YlOrRd") %>% colorRampPalette()
# NOTE that levelplot will depict the transposed transition matrix by default!
lattice::levelplot(t(ec_earth.CMIP5.transitionProb), col.regions = pcolors(201), at = seq(0,.6,0.05),
                  main = "Transition probabilities",
                   xlab = "", ylab = "",
                   scales = list(alternating = 3))

# Hypothesis Test for the Difference of Two Population Proportions ---------

pvalues <- transitionProb.test(obs.grid = clusters.era.interim, gcm.grid = wts.EC.earth)

lattice::levelplot(t(pvalues), col.regions = pcolors(201), at = seq(0,1,0.01),
                                      main = "EC-Earth prop.test - p.values",
                                      xlab = "", ylab = "",
                                      scales = list(alternating = 3))

# Significantly different transition probs are marked:
pvalues[which(pvalues >= 0.05)] <- NA
pvalues[which(pvalues < 0.05)] <- 1

lattice::levelplot(t(pvalues), col.regions = "black",
                   main = "EC-Earth - Significantly different probabilities",
                   xlab = "", ylab = "", colorkey = FALSE,
                   scales = list(alternating = 3))


  ## for barplot: http://myrcodes.blogspot.com/2014/10/overlaying-graphs-from-different.html
