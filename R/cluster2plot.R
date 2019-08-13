#     cluster2plot.R Cluster analysis of grid data
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

#'@title From Grid of Clusters to Plot
#'@description Performs subsetting in variable and member dimensions on a grid of clusters and prepares it to be plotted.  
#'@param cluster A grid (gridded or station dataset), multigrid, multimember grid or multimember multigrid object, as 
#' returned e.g. by \code{clusterGrid} that contains clusters from a weather type in the time dimension.
#'@param members Integer value indicating \strong{the position} of the member to be subset. members=1 by default.
#'If input grid has no member dimension, this argument is ignored. 
#'@param var Character vector indicating the variables(s) to be extracted. (Used for multigrid subsetting). 
#'See \link[transformeR]{subsetGrid}. It takes the first variable by default. 
#'@seealso \link[transformeR]{makeMultiGrid}, \link[transformeR]{subsetGrid}, \link[transformeR]{subsetDimension}.
#'@return A new grid object that contains the clusters from the specified variable and/or member. 
#'This subset grid of clusters is ready to be plotted with Climate4R plotting tools, e.g. \code{spatialPlot}.
#'@details  
#'@author J. A. Fernandez
#'@export
#'@importFrom transformeR makeMultiGrid subsetGrid subsetDimension getVarNames
#'@examples #Example 1: 'cluster' is a 3D grid of clusters.
#'cluster<- clusterGrid(SLP, type="kmeans", centers=10, iter.max=1000)
#'mg <- cluster2plot(cluster)
#'spatialPlot(mg, backdrop.theme = "coastline", rev.colors = T, layout = c(2,ceiling(attr(clusters, "centers")/2)), as.table=TRUE)
#'
#'Example 2: 'cluster' is a grid of clusters with variables and/or members. 
#'clusters<- clusterGrid(SLP.complex.var, type="kmeans", centers=10, iter.max=1000)
#'mg <- cluster2plot(clusters, members=1, var="slp") #To obtain first member and variable called "slp"
#'spatialPlot(mg, backdrop.theme = "coastline", rev.colors = T, layout = c(2,ceiling(attr(clusters, "centers")/2)), as.table=TRUE)



cluster2plot <- function(cluster, members=1, var=getVarNames(cluster)[1]){
  
  #check if input grid has cluster_type attrib. 
  if(is.null(attr(cluster, "cluster.type"))){
    stop("Input grid is not a grid of clusters.")
  }
  s<- subsetGrid(cluster, var=var, members = members, drop=TRUE) #dimension members was = 2, after this is = 1
  cluster.grids <- lapply(1:attr(s, "centers"), function(x) {
    subsetDimension(s, dimension="time", indices=x)})
  mg <- do.call("makeMultiGrid", c(cluster.grids, skip.temporal.check = TRUE))
  attr(mg$Variable,"longname")<-paste("Cluster_", 1:getShape(mg, "var"),sep="") #Prints Cluster_1, Cluster_2, etc...
  return(mg)
}

