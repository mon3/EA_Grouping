library(cec2005benchmark)
library(cec2013)
library(DEoptim)
library(genalg)
library(GA)
library(cluster)
library(clValid)
library(parallel)
library(doParallel)
library(dbscan)


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

pamClustering <- function(dataset, minimum=2, maximum){
  
  distance <- dist(dataset, method = "euclidean")
  bestGrouping = as.data.frame(rep(1, 100))
  bestWidth = -1.0
  
  for(i in minimum:maximum){
    fit <- pam(dataset, i, diss = inherits(dataset, "dist"), metric = "euclidean")
    currentGrouping <- fit$clustering
    currentSil <- fit$silinfo
    currentWidth <- fit$silinfo$avg.width
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



# finds optimal eps parameter for the dataset depending on the allowed % of noise points
# returns vector of clusters with optimal eps parameter
dbscanAnalyse <- function(dataset, minPts, dim=100){
  part = 0.5;
  i = 1.0
  while (part > 0.2)
  {
   # print("here")
    ds <- dbscan(dataset, eps = i, minPts = 11)
    i = i+1.0
    part = sum(ds$cluster==0)/100
  }
  print (ds$eps)
  return (ds$cluster)
}

GA.pop = c()
DE.pop = c()
RBGA.pop = c()
results = list()

for(funNr in 7:9){
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
  resultList = list(GA = list(), DE = list(), RBGA = list())
  
  for(i in 0:9){
    start = i*100+1
    end = i*100+100
    GA.current <- GA.cldata[c(start:end),]
    DE.current <- DE.cldata[c(start:end),]
    RBGA.current <- RBGA.cldata[c(start:end),]
    GA.dist <- dist(GA.current)
    DE.dist <- dist(DE.current)
    RBGA.dist <- dist(RBGA.current)
    
    GAres = list()
    DEres = list()
    RBGAres = list()
    GA.Hgroup <- getBestHClust(2, 15, GA.current, "euclidean")
    # jako ze mamy wektory cech 10-wym, to z wykresów za dużo nie wynika
    # plot(GA.current, GA.Hgroup)
    GA.Kgroup <- kMeansClustering(GA.current, 2, 15)
    GA.Pgroup <- pamClustering(GA.current, 2, 15)
    #GA.Dgroup <- dbscanAnalyse(GA.current, 11, 100)
    
    DE.Hgroup <- getBestHClust(2, 15, DE.current, "euclidean")
    DE.Kgroup <- kMeansClustering(DE.current, 2, 15)
    DE.Pgroup <- pamClustering(DE.current, 2, 15)
    #DE.Dgroup <- dbscanAnalyse(DE.current, 11, 100)
    
    RBGA.Hgroup <- getBestHClust(2, 15, RBGA.current, "euclidean")
    RBGA.Kgroup <- kMeansClustering(RBGA.current, 2, 15)
    RBGA.Pgroup <- pamClustering(RBGA.current, 2, 15)
    #RBGA.Dgroup <- dbscanAnalyse(RBGA.current, 11, 100)
    
    
    
    GAres$Kdunn <- dunn(GA.dist, GA.Hgroup)
    GAres$Pdunn <- dunn(GA.dist, GA.Pgroup)
    GAres$Hdunn <- dunn(GA.dist, GA.Kgroup)
  #  GAres$Ddunn <- dunn(GA.dist, GA.Dgroup)
    DEres$Kdunn <- dunn(DE.dist, DE.Hgroup)
    DEres$Pdunn <- dunn(DE.dist, DE.Pgroup)
    DEres$Hdunn <- dunn(DE.dist, DE.Kgroup)
   # DEres$Ddunn <- dunn(DE.dist, DE.Dgroup)
    RBGAres$Kdunn <- dunn(RBGA.dist, RBGA.Hgroup)
    RBGAres$Pdunn <- dunn(RBGA.dist, RBGA.Pgroup)
    RBGAres$Hdunn <- dunn(RBGA.dist, RBGA.Kgroup)
  #  RBGAres$Ddunn <- dunn(RBGA.dist, RBGA.Dgroup)
    GAres$Ksil <- summary(silhouette(GA.Kgroup, GA.dist))$avg.width
    GAres$Psil <- summary(silhouette(GA.Pgroup, GA.dist))$avg.width # mozna by wyciagac silhouette z samego pam, ale nalezaloby przerobic wiecej rzeczy przed tym krokiem
    GAres$Hsil <- summary(silhouette(GA.Hgroup, GA.dist))$avg.width
  #  GAres$Dsil <- summary(silhouette(GA.Dgroup, GA.dist))$avg.width
    DEres$Ksil <- summary(silhouette(DE.Kgroup, DE.dist))$avg.width
    DEres$Psil <- summary(silhouette(DE.Pgroup, DE.dist))$avg.width
    DEres$Hsil <- summary(silhouette(DE.Hgroup, DE.dist))$avg.width
   # DEres$Dsil <- summary(silhouette(DE.Dgroup, DE.dist))$avg.width
    RBGAres$Ksil <- summary(silhouette(RBGA.Kgroup, RBGA.dist))$avg.width
    RBGAres$Psil <- summary(silhouette(RBGA.Pgroup, RBGA.dist))$avg.width
    RBGAres$Hsil <- summary(silhouette(RBGA.Hgroup, RBGA.dist))$avg.width
   # RBGAres$Dsil <- summary(silhouette(RBGA.Dgroup, RBGA.dist))$avg.width
    maxKdunn = max(GAres$Kdunn, DEres$Kdunn, RBGAres$Kdunn)
    maxHdunn = max(GAres$Hdunn, DEres$Hdunn, RBGAres$Hdunn)
    maxPdunn = max(GAres$Pdunn, DEres$Pdunn, RBGAres$Pdunn)
   # maxDdunn = max(GAres$Ddunn, DEres$Ddunn, RBGAres$Ddunn)
    maxKsil = max(GAres$Ksil, DEres$Ksil, RBGAres$Ksil)
    maxHsil = max(GAres$Hsil, DEres$Hsil, RBGAres$Hsil)
    maxPsil = max(GAres$Psil, DEres$Psil, RBGAres$Psil)
  #  maxDsil = max(GAres$Dsil, DEres$Dsil, RBGAres$Dsil)
    resultList$GA[[i+1]] = GAres
    resultList$DE[[i+1]] = DEres
    resultList$RBGA[[i+1]] = RBGAres
    
    if(maxKdunn==GAres$Kdunn){
      GA.bestCounter <- GA.bestCounter + 1
    } else if(maxKdunn==DEres$Kdunn){
      DE.bestCounter <- DE.bestCounter + 1
    } else if(maxKdunn==RBGAres$Kdunn){
      RBGA.bestCounter <- RBGA.bestCounter + 1
    }
    
    if(maxHdunn==GAres$Hdunn){
      GA.bestCounter <- GA.bestCounter + 1
    } else if(maxHdunn==DEres$Hdunn){
      DE.bestCounter <- DE.bestCounter + 1
    }  else if(maxHdunn==RBGAres$Hdunn){
      RBGA.bestCounter <- RBGA.bestCounter + 1
    }
    
    if(maxPdunn==GAres$Pdunn){
      GA.bestCounter <- GA.bestCounter + 1
    } else if(maxPdunn==DEres$Pdunn){
      DE.bestCounter <- DE.bestCounter + 1
    }  else if(maxPdunn==RBGAres$Pdunn){
      RBGA.bestCounter <- RBGA.bestCounter + 1
    }
    
    # if(maxDdunn==GA.Ddunn){
    #   GA.bestCounter <- GA.bestCounter + 1
    # } else if(maxDdunn==DE.Ddunn){
    #   DE.bestCounter <- DE.bestCounter + 1
    # }  else if(maxDdunn==RBGA.Ddunn){
    #   RBGA.bestCounter <- RBGA.bestCounter + 1
    # }
    # 
    if(maxKsil==GAres$Ksil){
      GA.bestCounter <- GA.bestCounter + 1
    } else if(maxKsil==DEres$Ksil){
      DE.bestCounter <- DE.bestCounter + 1
    } else if(maxKsil==RBGAres$Ksil){
      RBGA.bestCounter <- RBGA.bestCounter + 1
    }
    
    if(maxHsil==GA.Hsil){
      GA.bestCounter <- GA.bestCounter + 1
    } else if(maxHsil==DE.Hsil){
      DE.bestCounter <- DE.bestCounter + 1
    } else if(maxHsil==RBGA.Hsil){
      RBGA.bestCounter <- RBGA.bestCounter + 1
    }
    
    if(maxPsil==GAres$Psil){
      GA.bestCounter <- GA.bestCounter + 1
    } else if(maxPsil==DEres$Psil){
      DE.bestCounter <- DE.bestCounter + 1
    } else if(maxPsil==RBGAres$Psil){
      RBGA.bestCounter <- RBGA.bestCounter + 1
    }
    
    # if(maxDsil==GA.Dsil){
    #   GA.bestCounter <- GA.bestCounter + 1
    # } else if(maxDsil==DE.Dsil){
    #   DE.bestCounter <- DE.bestCounter + 1
    # } else if(maxDsil==RBGA.Dsil){
    #   RBGA.bestCounter <- RBGA.bestCounter + 1
    # }
    # 
  }
  
  maxBest = max(GA.bestCounter, DE.bestCounter, RBGA.bestCounter)
  
  if(maxBest==GA.bestCounter){
    print("GA")
  } else if(maxBest==DE.bestCounter){
    print("DE")
  } else if(maxBest==RBGA.bestCounter){
    print("RBGA")
  }
  results[[funNr]] = resultList
}
