#'
#'@title Cluster Grid Implementation
#'@description This function creates a grid of clusters from another grid. The user will choose the clustering algorithm.  
#'@param grid The input grid or station data to be subset. This is either a grid (or station data), as 
#' returned e.g. by \code{loadeR::loadGridData} (or \code{loadeR::loadStationData}), a
#' multigrid, as returned by \code{makeMultiGrid}, or other types of multimember grids
#' (possibly multimember grids) as returned e.g. by \code{loadeR.ECOMS::loadECOMS}.
#'@param type Selects the clustering algorithm between k-means and hierarchical clustering. 
#' The possible values for 'type' are "\strong{kmeans}", "\strong{hierarchical}", "\strong{SOM}", or NULL. It chooses
#' K-means by default. See 'Details' for further information.
#'@param centers The number of clusters, \strong{k}, or center points. If SOM clustering is choosen, a two-component 
#'vector of integers will be needed. See Details. 
#'@param iter.max the maximum number of iterations allowed for K-means algorithm 
#'@param nstart (for K-means algorithm) if centers is a number, how many random sets should be chosen?
#'@param method the agglomeration method to be used in hierarchical. This should be one of "ward.D", "ward.D2", "single",
#'"complete", "average", "mcquitty", "median" or "centroid". It chooses "complete" by default. 
#'@seealso Clustering Algorithm Help: \link[stats]{kmeans}, \link[stats]{hclust}, \link[kohonen]{som}.
#'@return A new grid object that contains the clusters created using the specified algorithm. 
#'@details While using Hierarchical (check \link[stats]{hclust} for further information) clusterGrid() allows the 
#'user either to especified the numbers of clusters to be obtained or not. In the case that the user does not set 
#'input 'centers' to a integer value, that is to determine the number of clusters, this will be calculated  
#'automatically taking into account quantiles criteria. The function \link[stats]{quantile} is used to obtain 
#'quartiles at 0.25 and 0.75 and this is used as the height-difference threshold where to cut the tree into clusters 
#'by using \link[stats]{cutree}. If the user sets the numbers of clusters, clusterGrid() prioritize user's decision. 
#'While using SOM (check \link[kohonen]{som} for further information) the function calculates 
#'48 clusters (8x6) by default with rectangular topology. The user can modified the number of clusters obtained
#'by passing a two-component vector of integers to the input argument 'centers'. 
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
    kmModel <- kmeans(grid.2D, centers, iter.max = iter.max, nstart =nstart) #Datos de entrenamiento en KNN     
    #Going back from 2D to 3D:
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
    if(length(centers)!=2){
      stop("in 'centers'.\n Unexpected lenght for this argument while using SOM. It must be a vector of 2 elements")
    }else{
      som.grid<-som(grid.2D, somgrid(xdim=centers[1],ydim=centers[2],topo="rectangular"))
      Y <- mat2Dto3Darray(som.grid$codes[[1]], grid$xyCoords$x, grid$xyCoords$y)
    }
    
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