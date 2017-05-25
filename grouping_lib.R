library(cec2005benchmark)
library(cec2013)
library(DEoptim)
library(genalg)
library(GA)
library(cluster)
library(clValid)
library(parallel)
library(doParallel)

partial <- function(f, ...) {
  l <- list(...)
  function(...) {
    do.call(f, c(l, list(...)))
  }
}

gaSavePopulation <- function(obj, name){
  env = globalenv()
  env[[name]] = append(env[[name]], list(obj@population))
}

rbgaSavePopulation <- function(obj, name){
  env = globalenv()
  env[[name]] = append(env[[name]], list(obj$population))
}

populationToData <- function(obj){
  return(do.call(rbind.data.frame, obj))
}

#Zwraca wynik standardowego hclusta z podzia³em na k klastrów
#w zakresie minimum:maximum, gdzie k jest wybierane na podstawie
#najwy¿szej wartoœci œredniej silhouette dla zbioru
#metric - dowolna metryka, któr¹ przyjmuje dist
#lub "squared" dla kwadratowej euklidesowej
getBestHClust <- function(minimum = 2, maximum, dataset, metric){
  if(metric=="squared"){
    distance <- dist(dataset, method = "euclidean")
    distance <- distance^2
  }
  else
    distance <- dist(dataset, method = metric)
  
  fit = agnes(distance)
  bestGrouping = as.data.frame(rep(1, 100))
  bestWidth = -1.0
  for(i in minimum:maximum){
    currentGrouping <- cutree(fit, k=i)
    currentSil <- silhouette(currentGrouping, distance)
    currentWidth <- summary(currentSil)$avg.width
    if(currentWidth>bestWidth){
      bestGrouping <- currentGrouping
      bestWidth <- currentWidth
    }
  }
  return(bestGrouping)
}

#Zwraca odsetek grup, dla których wartoœæ funkcji dla œredniej jest lepsza (ni¿sza)
#ni¿ minimalna wartoœæ funkcji dla elementów grupy
groupIndex <- function(dataset, grouping, fun){
  subsets = c(c(0))
  groupsNr = max(grouping)
  betterNr = 0
  for(i in 1:groupsNr){
    currentSubset = data.matrix(subset(dataset, grouping==i))
    subsetMean = colMeans(currentSubset)
    if(min(fun(currentSubset))>fun(subsetMean)){
      betterNr <- betterNr + 1
    }
  }
  return(betterNr/groupsNr)
}