#     clusterGrid.R Cluster analysis of grid data
#
#     Copyright (C) 2019 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#'@title Cluster analysis of grids
#'@description Performs cluster analysis of grids, multigrids or multimember multigrids. Several clustering algorithms are available.  
#'@param grid A grid (gridded or station dataset), multigrid, multimember grid or multimember multigrid object, as 
#' returned e.g. by \code{loadeR::loadGridData} (or \code{loadeR::loadStationData}), a
#' multigrid, as returned by \code{makeMultiGrid}, or other types of multimember grids
#' (possibly multimember grids) as returned e.g. by \code{loadeR.ECOMS::loadECOMS}.
#'@param type Clustering algorithm to be used for the cluster analysis. 
#'Possible values are "\strong{kmeans}" (default), "\strong{hierarchical}", "\strong{som}". 
#'The core functions are \link[stats]{kmeans}, \link[stats]{hclust}, \link[kohonen]{som}, respectively. See Details.
#'@param centers Integer value indicating the number of clusters, \strong{k}, or center points. See Details.
#'@param iter.max (for the K-means algorithm) Integer value indicating the maximum number of iterations allowed. Default: 10.
#'@param nstart (for the K-means algorithm) If centers is a number, how many random sets should be chosen? Default: 1.
#'@param method (for the hierarchical algorithm) Agglomeration method to be used, one of "complete" (default), "ward.D", "ward.D2", "single",
#'"average", "mcquitty", "median" or "centroid". 
#'@seealso \link[stats]{kmeans}, \link[stats]{hclust}, \link[kohonen]{som}.
#'@return A new grid object that contains the clusters created using the specified algorithm.
#'The clustering type, number of clusters and other algorithm-specific parameters are provided as attributes.
#'@details 
#'\strong{kmeans}
#'
#'While using the  K-means algorithm, the number of clusters (argument 'centers') needs to be provided (no default). 
#'See  \link[stats]{kmeans} for more details in the implementation.
#'
#'\strong{hierarchical}
#'
#'While using the hierarchical algorithm (check \link[stats]{hclust} for further information) 
#'\code{clusterGrid} allows the user either to especify the number of clusters ('centers') or not. 
#'If the argument 'centers' is not provided, they are automatically set and the tree is cut when the height 
#'difference between two consecutive divisions (sorted in ascending order) is larger than the intercuartile 
#'range of the heights vector (see \link[stats]{cutree}) . 
#'
#'\strong{som}
#'
#'While using the SOM (self-organizing maps) algorithm (check \link[kohonen]{som} for further information), the argument 'centers' is provided as
#' a two-element vector, indicating the dimensions \code{xdim,ydim} of the grid (see \link[kohonen]{somgrid}).
#'Otherwise, by default 48 clusters (8x6) with rectangular topology are obtained. 
#'@author J. A. Fernandez
#'@export
#'@importFrom transformeR array3Dto2Dmat mat2Dto3Darray
#'@importFrom stats kmeans hclust cutree
#'@importFrom kohonen som
#'@examples #Example of K-means clustering: 
#'grid<-loadGridData(dataset , var = "grid")
#'clusters<- clusterGrid(grid, type="kmeans", centers=10, iter.max=1000)
#'cluster.grids <- lapply(1:10, function(i) {
#'   subsetDimension(clusters, dimension="time", indices=i)})
#'mg <- do.call("makeMultiGrid", c(cluster.grids, skip.temporal.check = TRUE))
#'spatialPlot(mg, backdrop.theme = "coastline", rev.colors = T, layout = c(2,5))
#'#Example of hierarchical clustering: 
#'...
#'#Example of som clustering: 
#'...


clusterGrid <- function(grid, type="kmeans", centers=NULL, iter.max=10, nstart=1, method = "complete"){
  
  #Argumento Type: para distinguir entre Kmeans, jerarquico, som (redes neuronales  ). 
  #if type=empty -> Kmeans by default
  
  type=tolower(type)
  grid.2D <- array3Dto2Dmat(grid$Data) #From 3D to 2D. pasamos a 2D ya qyue Kmeans trabaja con matrices.
  
  if(is.null(type) | type == "kmeans"){
    kmModel <- kmeans(grid.2D, centers, iter.max = iter.max, nstart =nstart) #Datos de entrenamiento en KNN     
    #Going back from 2D to 3D:
    Y <- mat2Dto3Darray(kmModel$centers, grid$xyCoords$x, grid$xyCoords$y)
    
  }else if (type == "hierarchical"){
    hc <- hclust(dist(grid.2D), method)
    if(is.null(centers)){
      #Auto-calculation of the number of clusters: Quartile method applied
      quantile.range<-quantile(hc$height, c(0.25,0.75))
      hc.height.diff<-numeric(length(hc$height)-1)
      for (i in 1:length(hc$height)-1){
        hc.height.diff[i]<-hc$height[i+1]-hc$height[i]
      }
      index <- which(hc.height.diff > (quantile.range[[2]]-quantile.range[[1]]))
      centers <- length(hc$order)-index[[1]]
    }
    memb <- cutree(hc, k = centers) #Found the corresponding cluster for every element
    cent <- NULL
    for(k in 1:centers){   #Set up the centers of the clusters
      cent <- rbind(cent, colMeans(grid.2D[memb == k, , drop = FALSE]))
    }
    Y <- mat2Dto3Darray(cent, grid$xyCoords$x, grid$xyCoords$y)
    
  }else if (type == "som"){
    if(is.null(centers)){
      som.grid<-som(grid.2D)
    }else {
      if(length(centers)!=2){
        stop("in 'centers'.\n Unexpected lenght for this argument while using SOM. It must be a vector of 2 elements")
      }
      som.grid<-som(grid.2D, somgrid(xdim=centers[1],ydim=centers[2],topo="rectangular"))
    }
    Y <- mat2Dto3Darray(som.grid$codes[[1]], grid$xyCoords$x, grid$xyCoords$y)
  } else {
    stop("Input data is not valid.\n'", paste(type),"' is not a valid algorithm")
  }
  
  #Setting up metadata for Y
  aux <- grid
  aux$Data <- Y
  attr(aux, "cluster_type")<- type
  #Add heights for hierarchical as attribute
  if(type == "hierarchical"){
    attr(aux, "height")<- hc$height
  }
  
  aux$Variable$varName<-"clusters"

  return(aux)

}