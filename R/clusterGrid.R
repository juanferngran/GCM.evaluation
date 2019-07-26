#'
#'@title Cluster Grid Implementation
#'@description This function creates a grid of clusters from another grid. The user will choose the clustering algorithm.  
#'@param grid The input grid or station data to be subset. This is either a grid (or station data), as 
#' returned e.g. by \code{loadeR::loadGridData} (or \code{loadeR::loadStationData}), a
#' multigrid, as returned by \code{makeMultiGrid}, or other types of multimember grids
#' (possibly multimember grids) as returned e.g. by \code{loadeR.ECOMS::loadECOMS}.
#'@param type Selects the clustering algorithm between k-means and hierarchical clustering. 
#' The possible values for 'type' are "\strong{kmeans}", "\strong{hierarchical}", or NULL. It chooses K-means by default. 
#' See 'Details' for further information.
#'@param centers The number of clusters, \strong{k}, or center points.  
#'@param iter.max the maximum number of iterations allowed for K-means algorithm 
#'@param nstart (for K-means algorithm) if centers is a number, how many random sets should be chosen?
#'@seealso Clustering Algorthim Help: \link[stats]{kmeans} 
#'@return A new grid object that contains the clusters created using the specified algorithm. 
#'@details For further information see \link{kmeans} function.
#'@keywords 
#'@author J. A. Fernandez
#'@export
#'@importFrom transformeR getShape aggregateGrid isRegular
#'@examples #Example of the implementation of K-means clustering and their representation: 
#'grid<-loadGridData(dataset , var = "grid")
#'clusters<- clusterGrid(grid, kmeans, 10, 1000, 1)
#'cluster.grids <- lapply(1:10, function(i) {
#'   subsetDimension(clusters, dimension="time", indices=i)})
#'mg <- do.call("makeMultiGrid", c(cluster.grids, skip.temporal.check = TRUE))
#'spatialPlot(mg, backdrop.theme = "coastline", rev.colors = T, layout = c(2,5))


clusterGrid <- function(grid, type="kmeans", centers, iter.max=10, nstart=1, method = "complete"){
  
  #Argumento Type: para distinguir entre Kmeans, jerarquico, som (redes neuronales  ). 
  #if type=empty -> Kmeans by default
  
  type=tolower(type)
  grid.2D <- array3Dto2Dmat(grid$Data) #From 3D to 2D. pasamos a 2D ya qyue Kmeans trabaja con matrices.
  
  if(is.null(type) | type == "kmeans"){
    #Realizamos clusters en la dimensión tiempo:
    kmModel <- kmeans(grid.2D, centers, iter.max = iter.max, nstart =nstart) #Datos de entrenamiento en KNN      #Going back from 2D to 3D:
    Y <- mat2Dto3Darray(kmModel$centers, grid$xyCoords$x, grid$xyCoords$y)
    
  }else if (type == "hierarchical"){
    hc <- hclust(dist(grid.2D), method)
    if(is.null(centers)){
      #Auto-calculation of the number of clusters: Quartile method applied
      quantile.range<-quantile(hc$height, c(0.25,0.75))
      hc.height.diff<-numeric(length(hc$height)-1)
      for (i in 1:length(hc$height)){
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