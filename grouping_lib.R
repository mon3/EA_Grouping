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

#Zwraca wynik standardowego hclusta z podzia?em na k klastr?w
#w zakresie minimum:maximum, gdzie k jest wybierane na podstawie
#najwy?szej warto?ci ?redniej silhouette dla zbioru
#metric - dowolna metryka, kt?r? przyjmuje dist
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

#Zwraca odsetek grup, dla kt?rych warto?? funkcji dla ?redniej jest lepsza (ni?sza)
#ni? minimalna warto?? funkcji dla element?w grupy
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

# do wyznaczania optymalej liczby klastrów dla k-means, maxClusters - maksymalna dopuszczalna liczba klastrów
# obecnie: nie jest używana
determineNumberOfClusters <- function(dataset, maxClusters){
  wss <- (nrow(dataset)-1)*sum(apply(dataset,2,var))
  for (i in 2:maxClusters) wss[i] <- sum(kmeans(dataset,centers=i)$withinss)
  plot(1:maxClusters, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}

# tylko dla metryki Euklidesowej!
# wyznaczenie optymalnej liczby klastrów dla grupowania k-średnich
# na podstawie współczynnika Silhouette
kMeansClustering <- function(dataset, minimum=2, maximum){
  
  distance <- dist(dataset, method = "euclidean")
  bestGrouping = as.data.frame(rep(1, 100))
  bestWidth = -1.0
  
  for(i in minimum:maximum){
    fit <- kmeans(dataset, i)
    currentGrouping <- fit$cluster
    currentSil <- silhouette(currentGrouping, distance)
    currentWidth <- summary(currentSil)$avg.width
    if(currentWidth>bestWidth){
      bestGrouping <- currentGrouping
      bestWidth <- currentWidth
    }
  }
  aggregate(dataset, by = list(bestGrouping), FUN = mean)
  dataset <- data.frame(dataset, bestGrouping)
  return(bestGrouping)
}


# dla populacji 100 osobników, wybieramy co 1/10 iteracji alg. ewolucyjnego
populationToClusterAnalysis <- function(data){
  popToClust = matrix(0, nrow = 1000, ncol = 10)
  popToClust <- data.frame(popToClust)
  
  for (i in 1:10){
    ind1low = i*100-99
    ind1high=i*100
    ind2low = i*0.1*nrow(data)-99
    ind2high = i* 0.1 * nrow(data)
    popToClust[ind1low:ind1high,] = data[ind2low:ind2high,]
  }
  # zwraca zbiór populacji z poszczególnych kroków AE, dla której będzie przeprowadzona analiza grupowania
  return(popToClust)
}


