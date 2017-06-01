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
#distance - macierz odleglosci grupowanych danych
#type - metoda grupowania spo?r?d dopuszczalnych przez agnes
getBestHClust <- function(minimum = 2, maximum, distance, type = "average"){
  fit = agnes(distance, method = type)
  bestGrouping = as.data.frame(rep(1, 100))
  bestWidth = -1.0
  for(i in minimum:maximum){
    currentGrouping <- cutree(fit, k=i)
    currentSil <- silhouette(currentGrouping, distance)
    currentWidth <- summary(currentSil)$avg.width
    print("CURRENT WIDTH =")
    print(currentWidth)
    print("BEST WIDTH = ")
    print(bestWidth)
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
kMeansClustering <- function(distance, minimum=2, maximum){
  bestGrouping = as.data.frame(rep(1, 100))
  bestWidth = -1.0

  for(i in minimum:maximum){
    fit <- kmeans(distance, i, 20)
    currentGrouping <- fit$cluster
    currentSil <- silhouette(currentGrouping, distance)
    currentWidth <- summary(currentSil)$avg.width
    if(currentWidth>bestWidth){
      bestGrouping <- currentGrouping
      bestWidth <- currentWidth
    }
  }
  return(bestGrouping)
}

#tylko dla metryki euklidesowej lub Manhattan!
pamClustering <- function(distance, minimum=2, maximum){
  bestGrouping = as.data.frame(rep(1, 100))
  bestWidth = -1.0

  for(i in minimum:maximum){
    fit <- pam(distance, i)
    currentGrouping <- fit$clustering
    currentSil <- fit$silinfo
    currentWidth <- fit$silinfo$avg.width
    if(currentWidth>bestWidth){
      bestGrouping <- currentGrouping
      bestWidth <- currentWidth
    }
  }
  return(bestGrouping)
}





# dla populacji 100 osobnikow, wybieramy co 1/10 iteracji alg. ewolucyjnego
populationToClusterAnalysis <- function(data, popNr){
  popToClust = matrix(0, nrow = popNr*100, ncol = 10)
  popToClust <- data.frame(popToClust)
  coeff = 1/popNr
  print("ROWS = ")
  print(nrow(popToClust))
  print("Coeff= ")
  print(coeff)
  for (i in 1:popNr){
    ind1low = i*100-99
    ind1high=i*100
    ind2low = i*coeff*nrow(data)-99
    ind2high = i* coeff * nrow(data)
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

gaPopulation <- function(funNr, poplSize, indivSize, iterations, crossProb, mutProb){
  ga(type = "real-valued", fitness = partial(cec2013, i=funNr), min = rep(poplSize*-1, indivSize), max = rep(poplSize, indivSize),
     maxiter = iterations, popSize=poplSize, pcrossover = crossProb, pmutation = mutProb, parallel = TRUE, monitor = partial(gaSavePopulation, name="GA.pop.anal"), seed=12345)
}


# algType: PAM, k-means; indexType: 1-silhouette, 2-Dunn, 3-averages;
kAlgorithmGroupingQuality <- function(minK=2, maxK, algType, indexType){
  GA.anal.cldata <- populationToClusterAnalysis(populationToData(GA.pop.anal), 10) # populacje z poszczególnych kroków GA

  distance <- dist(GA.anal.cldata[901:1000,], method = "euclidean")
  bestGrouping = as.data.frame(rep(1, 100))
  bestWidth = -1.0
  resSil <- vector()
  resDun <- vector()

  if (algType == "PAM"){
    for(i in minK:maxK){
      fit <- pam(dataset, i, diss = inherits(dataset, "dist"), metric = "euclidean")
      currentGrouping <- fit$clustering
      currentSil <- fit$silinfo
      currentWidth <- fit$silinfo$avg.width
      resSil[i-1] <- fit$silinfo$avg.width
      resDun[i-1] <- dunn(distance, currentGrouping)
      # zapisujemy silhouette i k
      if(currentWidth>bestWidth){
        bestGrouping <- currentGrouping
        bestWidth <- currentWidth
      }
    }
  }
  else if (algType == "k-means"){

    for(i in minK:maxK){
      fit <- kmeans(dataset, i)
      currentGrouping <- fit$cluster
      currentSil <- silhouette(currentGrouping, distance)
      currentWidth <- summary(currentSil)$avg.width
      resSil[i-1] <- currentWidth
      resDun[i-1] <- dunn(distance, currentGrouping)
      # zapisujemy silhouette i k
      if(currentWidth>bestWidth){
        bestGrouping <- currentGrouping
        bestWidth <- currentWidth
      }
    }
  }

  if (indexType==1){
    print("SILHOUETTE")
    print (resSil)
    print("Optimal k = ")
    print(which(resSil==max(resSil)))

  }
  else if (indexType==2){
    print("DUNN")
    print(resDun)
    print("Optimal k = ")
    print(which(resDun==max(resDun)))
  }

}

assessGroupingAlgorithm <- function(data, npops, algorithm, func, hmethod = "average", metric = "euclidean", kMin = 2, kMax = 15, tryAllK = FALSE){
  groupingResults = list()
  step = nrow(data)/npops

  for(i in 0:(npops-1)){
    start = i*step+1
    end = i*step+100
    pop.current <- data[c(start:end),]
    if(algorithm == "hclust" && !is.null(metric)){
      if(metric == "squared"){
          pop.dist <- dist(pop.current, method = "euclidean")
          pop.dist <- pop.dist^2
      } else{
        pop.dist <- dist(pop.current, method = metric)
      }
    } else{
      pop.dist <- dist(pop.current)
    }

    partialRes = list()

    if(tryAllK == TRUE){
      for(k in kMin:kMax){
        kRes = list()
        if(algorithm == "hclust"){
          kRes$grouping = getBestHClust(k, k, pop.dist, hmethod)
        } else if(algorithm == "kclust"){
          kRes$grouping = kMeansClustering(pop.dist, k, k)
        } else if(algorithm == "pam"){
          kRes$grouping = pamClustering(pop.dist, k, k)
        }
        kRes$dunn = dunn(pop.dist, kRes$grouping)
        kRes$sil = summary(silhouette(kRes$grouping, pop.dist))$avg.width
        kRes$gind = groupIndex(pop.current, kRes$grouping, func)
        partialRes[[k]] = kRes
      }
    } else{
      kRes = list()
      if(algorithm == "hclust"){
        kRes$grouping = getBestHClust(kMin, kMax, pop.dist, hmethod)
      } else if(algorithm == "kclust"){
        kRes$grouping = kMeansClustering(pop.dist, kMin, kMax)
      } else if(algorithm == "pam"){
        kRes$grouping = pamClustering(pop.dist, kMin, kMax)
      }
      kRes$dunn = dunn(pop.dist, kRes$grouping)
      kRes$sil = summary(silhouette(kRes$grouping, pop.dist))$avg.width
      kRes$gind = groupIndex(pop.current, kRes$grouping, func)
      partialRes[[max(kRes$grouping)]] = kRes
    }
    groupingResults[[i+1]] = partialRes
  }
  return(groupingResults)
}

# returns list of results for different grouping algorithms
runGAGrouping <- function(GaGroupsNr, gaGroups, popNr, metrics, funNr, iterations, crossProb, mutProb, cecNr){
  GA.pop.anal = c()
  gaPopulation(funNr, 100, 10, iterations, crossProb, mutProb)
  GA.pop.anal.cldata <- populationToClusterAnalysis(populationToData(GA.pop.anal), popNr)
  
  fun = partial(cec2013, i=cecNr)
  
  groupingResPAM <- list()
  groupingResH <- list()
  groupingResK <- list()
  
  # GRUPOWANIE
  if (GaGroupsNr==3){
    groupingResK = assessGroupingAlgorithm(GA.pop.anal.cldata, popNr, "kclust", fun, tryAllK = TRUE)
    groupingResH = assessGroupingAlgorithm(GA.pop.anal.cldata, popNr, metric = metrics, "hclust", fun)
    groupingResPAM = assessGroupingAlgorithm(GA.pop.anal.cldata, popNr, "pam", fun, tryAllK = TRUE)
  } else if (GaGroupsNr==2){
    if (gaGroups[1]=="hclust" || gaGroups[2]== "hclust"){
      groupingResH = assessGroupingAlgorithm(GA.pop.anal.cldata, popNr, metric = metrics, "hclust", fun)
    }
    if (gaGroups[1]=="pam" || gaGroups[2]== "pam"){
      groupingResPAM = assessGroupingAlgorithm(GA.pop.anal.cldata, popNr, "pam", fun, tryAllK = TRUE)
    }
    if (gaGroups[1]=="kmeans" || gaGroups[2]== "kmeans"){
      groupingResK = assessGroupingAlgorithm(GA.pop.anal.cldata, popNr, "kclust", fun, tryAllK = TRUE)
    }
  } else
  {
    if (gaGroups[1]=="hclust" ){
      groupingResH = assessGroupingAlgorithm(GA.pop.anal.cldata, popNr, metric = metrics, "hclust", fun)
    }
    if (gaGroups[1]=="pam"){
      groupingResPAM = assessGroupingAlgorithm(GA.pop.anal.cldata, popNr, "pam", fun, tryAllK = TRUE)
    }
    if (gaGroups[1]=="kmeans"){
      groupingResK = assessGroupingAlgorithm(GA.pop.anal.cldata, popNr, "kclust", fun, tryAllK = TRUE)
    }
  }
  grouping <- c(groupingResK, groupingResH, groupingResPAM)
  return(grouping)
}



source("read_params.R")
parameters = readInputParams()
metrics = parameters$metrics
popNr = parameters$popNr
agnesType = parameters$agnesType
GaGroupsNr = parameters$GaGroupsNr
gaGroups = parameters$GaGroups
# liczba opcjonalnych parametrów + wartości
GAparamsNr = parameters$GAparamsNr
GAparams = parameters$GAparams
# te parametry do wywołania GA ewentualnie
if (GAparamsNr != 0){
  GAfunNr = strtoi(GAparams[1])
  GAiteratins = strtoi(GAparams[2])
  GAcrossover = as.numeric(GAparams[3])
  GAmutation = as.numeric(GAparams[4])
}
#print(GAparams[2])
print("AGNES TYPE = ")
print(agnesType)


GA.pop = c()
GA.pop.anal = c()
DE.pop = c()
RBGA.pop = c()
results = list()

for(funNr in 7){
  ga(type = "real-valued", fitness = partial(cec2013, i=funNr), min = rep(-100, 10), max = rep(100, 10),
     maxiter = 1000, popSize=100, parallel = TRUE, monitor = partial(gaSavePopulation, name="GA.pop"))
  rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=100, iters = 1000,
       mutationChance=NA, elitism=NA, monitorFunc=partial(rbgaSavePopulation, name="RBGA.pop"),
       evalFunc = partial(cec2013, i=funNr))
  registerDoSEQ()
  DE.pop <- DEoptim(partial(cec2013, i=funNr), rep(-100, 10), rep(100, 10),
                       DEoptim.control(storepopfrom = 0, trace=FALSE, parallelType=2, itermax=1000))$member$storepop
  GA.cldata <- populationToClusterAnalysis(populationToData(GA.pop), 10)
  DE.cldata <- populationToClusterAnalysis(populationToData(DE.pop), 10)
  RBGA.cldata <- populationToClusterAnalysis(populationToData(RBGA.pop), 10)

  GA.bestCounter <- 0
  DE.bestCounter <- 0
  RBGA.bestCounter <- 0
  resultList = list(GA = list(), DE = list(), RBGA = list())

  for(i in 0:(popNr-1)){
    start = i*100+1
    end = i*100+100
    GA.current <- GA.cldata[c(start:end),]
    print(nrow(GA.current))
    # TODO: tutaj jest cos zle :(
    DE.current <- DE.cldata[c(start:end),]
    RBGA.current <- RBGA.cldata[c(start:end),]
    GA.dist <- dist(GA.current)
    print(GA.dist)
    DE.dist <- dist(DE.current)
    RBGA.dist <- dist(RBGA.current)

    GAres = list()
    DEres = list()
    RBGAres = list()
    
  #  print("METRICS  = ")
 #   print(metrics)

    GA.Hgroup <- getBestHClust(2, 15, GA.current, agnesType)

    # jako ze mamy wektory cech 10-wym, to z wykresów za dużo nie wynika
    # plot(GA.current, GA.Hgroup)
    GA.Kgroup <- kMeansClustering(GA.current, 2, 15)
    GA.Pgroup <- pamClustering(GA.current, 2, 15)
    #GA.Dgroup <- dbscanAnalyse(GA.current, 11, 100)
    
    DE.Hgroup <- getBestHClust(2, 15, DE.current, agnesType)
    DE.Kgroup <- kMeansClustering(DE.current, 2, 15)
    DE.Pgroup <- pamClustering(DE.current, 2, 15)
    #DE.Dgroup <- dbscanAnalyse(DE.current, 11, 100)
    
    RBGA.Hgroup <- getBestHClust(2, 15, RBGA.current, agnesType)
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
    if(maxHdunn==Gares$Hdunn){
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

    if(maxHsil==GAres$Hsil){
      GA.bestCounter <- GA.bestCounter + 1
    } else if(maxHsil==DEres$Hsil){
      DE.bestCounter <- DE.bestCounter + 1
    } else if(maxHsil==RBGAres$Hsil){
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

sums = c()
for(i in 1:3){
  partial = c()
  for(j in 7:9){
    sum = c()
    for(k in 1:6){
        sum[k] = Reduce("+",resjson[[i]][[j]][[k]])
    }
    partial[[j]] = sum
  }
  sums[[i]] = partial
}
