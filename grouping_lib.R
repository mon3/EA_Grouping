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

#Zwraca wynik standardowego hclusta z podzialem na k klastrow
#w zakresie minimum:maximum, gdzie k jest wybierane na podstawie
#najwyzszej wartosci sredniej silhouette dla zbioru
#metric - dowolna metryka, ktora przyjmuje dist
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

#Zwraca odsetek grup, dla ktorych wartosc funkcji dla sredniej jest lepsza (nizsza)
#niz minimalna wartosc funkcji dla elementow grupy
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

# do wyznaczania optymalej liczby klastrow dla k-means, maxClusters - maksymalna dopuszczalna liczba klastrow
# obecnie: nie jest uzywana
determineNumberOfClusters <- function(dataset, maxClusters){
  wss <- (nrow(dataset)-1)*sum(apply(dataset,2,var))
  for (i in 2:maxClusters) wss[i] <- sum(kmeans(dataset,centers=i)$withinss)
  plot(1:maxClusters, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}

# tylko dla metryki Euklidesowej!
# wyznaczenie optymalnej liczby klastrow dla grupowania k-srednich
# na podstawie wspolczynnika Silhouette
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


# dla populacji 100 osobnikow, wybieramy co 1/10 iteracji alg. ewolucyjnego
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
  # zwraca zbior populacji z poszczegolnych krokow AE, dla ktorej bedzie przeprowadzona analiza grupowania
  return(popToClust)
}

for(funNr in 7:11){
  ga(type = "real-valued", fitness = partial(cec2013, i=funNr), min = rep(-100, 10), max = rep(100, 10),
     maxiter = 1000, popSize=100, parallel = TRUE, monitor = partial(gaSavePopulation, name="GA.pop"))
  rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=100, iters = 1000, 
       mutationChance=NA, elitism=NA, monitorFunc=partial(rbgaSavePopulation, name="RBGA.pop"),
       evalFunc = partial(cec2013, i=funNr))
  registerDoSEQ()
  DE.pop <- DEoptim(partial(cec2013, i=funNr), rep(-100, 10), rep(100, 10),
                       DEoptim.control(storepopfrom = 0, trace=FALSE, parallelType=2, itermax=1000))$member$storepop
  GA.cldata <- populationToClusterAnalysis(populationToData(GA.pop))
  DE.cldata <- populationToClusterAnalysis(populationToData(DE.pop))
  RBGA.cldata <- populationToClusterAnalysis(populationToData(RBGA.pop))
  
  GA.bestCounter <- 0
  DE.bestCounter <- 0
  RBGA.bestCounter <- 0
  
  for(i in 0:9){
    GA.current <- GA.cldata[c(i*100+1:i*100+100),]
    DE.current <- DE.cldata[c(i*100+1:i*100+100),]
    RBGA.current <- RBGA.cldata[c(i*100+1:i*100+100),]
    GA.dist <- dist(GA.current)
    DE.dist <- dist(DE.current)
    RBGA.dist <- dist(RBGA.current)
    GA.Hgroup <- getBestHClust(2, 15, currentGA, "euclidean")
    GA.Kgroup <- kMeansClustering(currentGA, 2, 15)
    DE.Hgroup <- getBestHClust(2, 15, currentDE, "euclidean")
    DE.Kgroup <- kMeansClustering(currentDE, 2, 15)
    RBGA.Hgroup <- getBestHClust(2, 15, currentRBGA, "euclidean")
    RBGA.Kgroup <- kMeansClustering(currentRBGA, 2, 15)
    GA.Kdunn <- dunn(GA.dist, GA.Hgroup)
    GA.Hdunn <- dunn(GA.dist, GA.Kgroup)
    DE.Kdunn <- dunn(DE.dist, DE.Hgroup)
    DE.Hdunn <- dunn(DE.dist, DE.Kgroup)
    RBGA.Kdunn <- dunn(RBGA.dist, RBGA.Hgroup)
    RBGA.Hdunn <- dunn(RBGA.dist, RBGA.Kgroup)
    GA.Ksil <- silhouette(GA.Kgroup, GA.dist)
    GA.Hsil <- silhouette(GA.Hgroup, GA.dist)
    DE.Ksil <- silhouette(DE.Kgroup, DE.dist)
    DE.Hsil <- silhouette(DE.Hgroup, DE.dist)
    RBGA.Ksil <- silhouette(RBGA.Kgroup, RBGA.dist)
    RBGA.Hsil <- silhouette(RBGA.Hgroup, RBGA.dist)

    maxKdunn = max(GA.Kdunn, DE.Kdunn, RBGA.Kdunn)
    maxHdunn = max(GA.Hdunn, DE.Hdunn, RBGA.Hdunn)
    maxKsil = max(GA.Ksil, DE.Ksil, RBGA.Ksil)
    maxHsil = max(GA.Hsil, DE.Hsil, RBGA.Hsil)
    
    if(maxKdunn==GA.Kdunn){
      GA.bestCounter <- GA.bestCounter + 1
    }
    else if(maxKdunn==DE.Kdunn){
      DE.bestCounter <- DE.bestCounter + 1
    }
    if(maxKdunn==RBGA.Kdunn){
      RBGA.bestCounter <- RBGA.bestCounter + 1
    }
    
    if(maxHdunn==GA.Hdunn){
      GA.bestCounter <- GA.bestCounter + 1
    }
    else if(maxHdunn==DE.Hdunn){
      DE.bestCounter <- DE.bestCounter + 1
    }
    if(maxHdunn==RBGA.Hdunn){
      RBGA.bestCounter <- RBGA.bestCounter + 1
    }
    
    if(maxKsil==GA.Ksil){
      GA.bestCounter <- GA.bestCounter + 1
    }
    else if(maxKsil==DE.Ksil){
      DE.bestCounter <- DE.bestCounter + 1
    }
    if(maxKsil==RBGA.Ksil){
      RBGA.bestCounter <- RBGA.bestCounter + 1
    }
    
    if(maxHsil==GA.Hsil){
      GA.bestCounter <- GA.bestCounter + 1
    }
    else if(maxHsil==DE.Hsil){
      DE.bestCounter <- DE.bestCounter + 1
    }
    if(maxHsil==RBGA.Hsil){
      RBGA.bestCounter <- RBGA.bestCounter + 1
    }
  }
  
  maxBest = max(GA.bestCounter, DE.bestCounter, RBGA.bestCounter)
  
  if(maxBest==GA.bestCounter){
    print("GA")
  }
  else if(maxBest==DE.bestCounter){
    print("DE")
  }
  if(maxBest==RBGA.bestCounter){
    print("RBGA")
  }
}
