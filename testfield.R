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
    #print("CURRENT WIDTH =")
    #print(currentWidth)
    #print("BEST WIDTH = ")
    #print(bestWidth)
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
    fit <- kmeans(distance, i, 50, nstart = 3)
    if(fit$ifault==4){
      fit <- kmeans(distance, i, 50, algorithm="MacQueen")
    }
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
  coeff = nrow(data)/popNr
  #print("ROWS = ")
  #print(nrow(popToClust))
  #print("Coeff= ")
  #print(coeff)
  for (i in 1:popNr){
    ind1low = (i-1)*100+1
    ind1high= ind1low + 99
    ind2low = i*coeff-99 # nie daje (i-1), zeby miec ostatnia populacje
    ind2high = ind2low + 99
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
        } else if(algorithm == "kmeans"){
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
      } else if(algorithm == "kmeans"){
        kRes$grouping = kMeansClustering(pop.dist, kMin, kMax)
      } else if(algorithm == "pam"){
        kRes$grouping = pamClustering(pop.dist, kMin, kMax)
      }
      kRes$dunn <- dunn(pop.dist, kRes$grouping)
      kRes$sil <- summary(silhouette(kRes$grouping, pop.dist))$avg.width
      kRes$gind <- groupIndex(pop.current, kRes$grouping, func)
      kRes$k <- max(kRes$grouping)
      partialRes <- kRes
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

for(funNr in 6:15){
  currName1 = paste("pop1.", funNr, sep="")
  currName2 = paste("pop2.", funNr, sep="")
  env = globalenv()
  env[[currName1]] <- list()
  env[[currName2]] <- list()
rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=100, iters = 1000,
     mutationChance=NA, elitism=NA, monitorFunc=partial(rbgaSavePopulation, name=currName1),
                                                        evalFunc = partial(cec2013, i=funNr))
rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=100, iters = 1000,
     mutationChance=NA, elitism=NA, monitorFunc=partial(rbgaSavePopulation, name=currName2),
                                                        evalFunc = partial(cec2013, i=funNr))
}

resultMaster <- list()
npops <- 50
hclustMetrics <- list("euclidean", "maximum", "manhattan")
pamMetrics <- list("euclidean", "manhattan")
for(funNr in 6:15){
  resultSingle <- list()
  currName1 = paste("pop1.", funNr, sep="")
  currName2 = paste("pop2.", funNr, sep="")
  env = globalenv()
  current_p1 <- env[[currName1]]
  current_p2 <- env[[currName2]]
  data_p1 <- populationToClusterAnalysis(populationToData(current_p1), npops)
  data_p2 <- populationToClusterAnalysis(populationToData(current_p2), npops)
  hclust <- list()
  for(m in 1:length(hclustMetrics)){
    hclustSingle <- list()
    hclustSingle[[1]] <- assessGroupingAlgorithm(data_p1, npops, "hclust", partial(cec2013, i=funNr), 
                                                 metric = hclustMetrics[[m]])
    hclustSingle[[2]] <- assessGroupingAlgorithm(data_p2, npops, "hclust", partial(cec2013, i=funNr), 
                                                 metric = hclustMetrics[[m]])
    hclust[[hclustMetrics[[m]]]] <- hclustSingle
  }
  resultSingle$hclust <- hclust
  
  pamc <- list()
  for(m in 1:length(pamMetrics)){
    pamcSingle <- list()
    pamcSingle[[1]] <- assessGroupingAlgorithm(data_p1, npops, "pam", partial(cec2013, i=funNr), 
                                                 metric = pamMetrics[[m]])
    pamcSingle[[2]] <- assessGroupingAlgorithm(data_p2, npops, "pam", partial(cec2013, i=funNr), 
                                                 metric = pamMetrics[[m]])
    pamc[[pamMetrics[[m]]]] <- pamcSingle
  }
  resultSingle$pam <- pamc
  
  kmean <- list()
  kmean[[1]] <- assessGroupingAlgorithm(data_p1, npops, "kmeans", partial(cec2013, i=funNr))
  kmean[[2]] <- assessGroupingAlgorithm(data_p2, npops, "kmeans", partial(cec2013, i=funNr))
  resultSingle$kmeans <- kmean
  resultMaster[[funNr]] <- resultSingle
}