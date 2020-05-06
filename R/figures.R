
# 1. Figure 1: BAR PLOT with LWTs comparison of the four reanalysis in all Seasons: ------------------------

















# 1. Transition probabilities -----------------------------------------------------

## Add NCEP
# ds <- "http://meteo.unican.es/tds5/dodsC/ncepReanalysis1/ncepReanalysis1_4xDaily.ncml"
# wmo <- 1981:2010
# lonLim = c(-45, 66)
# latLim = c(22, 73)
# 
# library(loadeR)
# source("~/workspace/jb")
# loginUDG(username, password)
# ncep <- loadGridData(dataset = ds,
#                      var = "slp",
#                      lonLim = lonLim,
#                      latLim = latLim,
#                      season = c(12,1:11),
#                      years = wmo,
#                      time = "DD", 
#                      aggr.d = "mean")
# clusters.ncep <- clusterGrid(ncep, type = "lamb")
# save(ncep, clusters.ncep, file = "juan_WORK/Data/clustering_ncep.RData")

# # Add ERA20
# ds <- "http://meteo.unican.es/tds5/dodsC/ecmwf/era20c_sfc_151.128_an.ncml"
# wmo <- 1981:2010
# lonLim = c(-45, 66)
# latLim = c(22, 73)
# 
# library(loadeR)
# source("~/workspace/jb")
# loginUDG(username, password)
# era20 <- loadGridData(dataset = ds,
#                      var = "Mean_sea_level_pressure_surface",
#                      lonLim = lonLim,
#                      latLim = latLim,
#                      season = c(12,1:11),
#                      years = wmo,
#                      time = "DD",
#                      aggr.d = "mean")
# clusters.era20 <- clusterGrid(era20, type = "lamb")
# save(era20, clusters.era20, file = "juan_WORK/Data/clustering_era20.RData")


## 1.1. Transition probability plots ------------------------------------------------

rm(list = ls())
source("scripts/tprobPlot.R")
source("scripts/transitionProb.R")
source("scripts/transitionProb.pvalue.R")

load("./juan_WORK/Data/clustering_eraInterim.RData", verbose = TRUE)
load("./juan_WORK/Data/clustering_jra.RData", verbose = TRUE)
load("./juan_WORK/Data/clustering_ncep.RData", verbose = TRUE)
load("./juan_WORK/Data/clustering_era20.RData", verbose = TRUE)

cairo_pdf("figs/reanalysis_transitions.pdf", width = 8.63, height = 9.06, onefile = TRUE)
tprobPlot(tprob.matrix = transitionProb(clusters.era.interim), title = "ERA-Interim")
tprobPlot(tprob.matrix = transitionProb(clusters.jra), title = "JRA")
tprobPlot(tprob.matrix = transitionProb(clusters.ncep), title = "NCEP")
tprobPlot(tprob.matrix = transitionProb(clusters.ncep), title = "ERA-20C")
dev.off()


## 1.2. Observational uncertainty ---------------------------------------------------
tprob.erain <- transitionProb(clusters.era.interim)
tprob.jra <- transitionProb(clusters.jra)
tprob.ncep <- transitionProb(clusters.ncep)
tprob.era20 <- transitionProb(clusters.era20)

wtnames <- colnames(tprob.erain)


cairo_pdf("figs/missing_reanalysis_transitions.pdf", width = 7.67, height = 9.06)

par(mar = c(9, 4.5, 9, 5))
image(t(tprob.erain), col = NULL, axes = FALSE, main = "Non-observed transitions")
axis(1, at = seq(0, 1, 1/25), labels = wtnames, las = 2)
axis(2, at = seq(0, 1, 1/25), labels = wtnames, las = 1)
axis(3, at = seq(0, 1, 1/25), labels = wtnames, las = 2)
axis(4, at = seq(0, 1, 1/25), labels = wtnames, las = 1)
# fields::image.plot(tprob.matrix, legend.only = TRUE, col = pcolors(21))
coords <- seq(0, 1, 1/25)
## Missing in ERA-Interim
missing.ei <- which(is.na(tprob.erain), arr.ind = TRUE)
for (i in 1:nrow(missing.ei)) {
    points(coords[missing.ei[i,2]], coords[missing.ei[i,1]],
           pch = 21, cex = 1.75)
}
## Missing in JRA
missing.jra <- which(is.na(tprob.jra), arr.ind = TRUE)
for (i in 1:nrow(missing.jra)) {
    points(coords[missing.jra[i,2]], coords[missing.jra[i,1]],
           pch = 21, cex = 1.25, col = "red")
}
## Missing in NCEP
missing.ncep <- which(is.na(tprob.ncep), arr.ind = TRUE)
for (i in 1:nrow(missing.ncep)) {
    points(coords[missing.ncep[i,2]], coords[missing.ncep[i,1]],
           pch = 19, cex = .75, col = "black")
}
## Missing in ERA-20C
missing.e20 <- which(is.na(tprob.era20), arr.ind = TRUE)
for (i in 1:nrow(missing.e20)) {
    points(coords[missing.e20[i,2]], coords[missing.e20[i,1]],
           pch = 3, cex = 1, col = "blue")
}

par(xpd = TRUE) # Allow putting legend outside graphic area
legend(0, -0.15, ncol = 2,
       c(paste0("Does not occur in ERA-Interim (n=", nrow(missing.ei),")"), 
         paste0("Does not occur in JRA (n=", nrow(missing.jra), ")"),
         paste0("Does not occur in NCEP (n=", nrow(missing.ncep), ")"),
         paste0("Does not occur in ERA-20C (n=", nrow(missing.e20), ")")),
       pch = c(21, 21, 19, 3),
       pt.cex = c(1.75, 1.25, 0.75, 1),
       col = c("black", "red", "black", "blue"),
       cex = 1, bty = "n")

dev.off()

nrow(missing.ei)
nrow(missing.ncep)
nrow(missing.jra)

ind.erain <- which(is.na(tprob.erain))
ind.ncep <- which(is.na(tprob.ncep))
ind.jra <- which(is.na(tprob.jra))



# 2. GCM transitions -----------------------------------------------------------
rm(list = ls())

load("./juan_WORK/Data/clustering_eraInterim.RData", verbose = TRUE)

source("scripts/tprobPlot.R")
source("scripts/transitionProb.R")
source("scripts/transitionProb.pvalue.R")

## 2.1. transitions CMIP5 -----------------------------------------------------------


load("./juan_WORK/Data/WTs_cmip5.RData", verbose = TRUE)


# probatina: 
# tprobPlot(tprob.matrix = transitionProb(wts.EC.earth),
#           pval.matrix = transitionProb.test(obs.grid = clusters.era.interim,
#                                             gcm.grid =  wts.EC.earth),
#           tprob.ref = transitionProb(clusters.era.interim), title = "CMIP5 EC-EARTH") %>% print()

models <- c("wts.cccma", "wts.cnrm", "wts.EC.earth",
            "wts.ipsl", "wts.miroc", "wts.mohc",
            "wts.mpi.esm.lr", "wts.mpi.esm.mr",
            "wts.ncc.nor", "wts.noaa.gfdl" )

modelnames <- c("CANESM2", "CNRM-CM5", "EC-EARTH", "IPSL-CM5A-MR", "MIROC-ESM",
                "HADGEM2-ES", "MPI-ESM-LR", "MPI-ESM-MR", "NORESM-M1", "NOAA-GFDL-ESM2M") 

cairo_pdf("figs/CMIP5_transitions.pdf", width = 8.63, height = 9.06, onefile = TRUE)
for (i in 1:length(models)) {
    model <- get(models[i])
    tprobPlot(tprob.matrix = transitionProb(model),
              pval.matrix = transitionProb.test(obs.grid = clusters.era.interim,
                                                gcm.grid = model),
              tprob.ref = transitionProb(clusters.era.interim),
              title = paste0("CMIP5_", modelnames[i]))# %>% print()
}
dev.off()

## 2.2. transitions CMIP6 -----------------------------------------------------------

load("./juan_WORK/Data/WTs_cmip6.RData", verbose = TRUE)
models <- ls(pattern = "^wts\\..*\\.cmip6$")

modelnames <- c("CANESM5", "CNRM-CM6-1", "IPSL-CM6A-LR", "MIROC6",
                "HADGEM3-GC31-LL", "MPI-ESM1-2-LR", "NORESM2-LM") 

cairo_pdf("figs/CMIP6_transitions.pdf", width = 8.63, height = 9.06, onefile = TRUE)
for (i in 1:length(models)) {
    model <- get(models[i])
    tprobPlot(tprob.matrix = transitionProb(model),
              pval.matrix = transitionProb.test(obs.grid = clusters.era.interim,
                                                gcm.grid = model),
              tprob.ref = transitionProb(clusters.era.interim),
              title = paste0("CMIP6_", modelnames[i]))# %>% print()
}
dev.off()



