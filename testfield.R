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

cutAverageResHclust <- function(resultMaster){
  hscores = list()
  
  for(i in 6:15){
    singlescore = list()
    singlescoreMax = list()
    singlescoreMan = list()
    for(j in seq(1,50)){
      singlescore[[j]] <- list()
      singlescore[[j]]$dunn <- 0
      singlescore[[j]]$sil <- 0
      singlescore[[j]]$gind <- 0
      singlescore[[j]]$k <- 0
      
      singlescoreMax[[j]] <- list()
      singlescoreMax[[j]]$dunn <- 0
      singlescoreMax[[j]]$sil <- 0
      singlescoreMax[[j]]$gind <- 0
      singlescoreMax[[j]]$k <- 0
      
      singlescoreMan[[j]] <- list()
      singlescoreMan[[j]]$dunn <- 0
      singlescoreMan[[j]]$sil <- 0
      singlescoreMan[[j]]$gind <- 0
      singlescoreMan[[j]]$k <- 0
      for(k in 1:10){
        singlescore[[j]]$dunn <- singlescore[[j]]$dunn + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$dunn
        singlescore[[j]]$sil <- singlescore[[j]]$sil + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$sil
        singlescore[[j]]$gind <- singlescore[[j]]$gind+ resultMaster[[i]]$hclust$euclidean[[k]][[j]]$gind
        singlescore[[j]]$k <- singlescore[[j]]$k + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$k
        
        singlescoreMax[[j]]$dunn <- singlescoreMax[[j]]$dunn + resultMaster[[i]]$hclust$maximum[[k]][[j]]$dunn
        singlescoreMax[[j]]$sil <- singlescoreMax[[j]]$sil + resultMaster[[i]]$hclust$maximum[[k]][[j]]$sil
        singlescoreMax[[j]]$gind <- singlescoreMax[[j]]$gind+ resultMaster[[i]]$hclust$maximum[[k]][[j]]$gind
        singlescoreMax[[j]]$k <- singlescoreMax[[j]]$k + resultMaster[[i]]$hclust$maximum[[k]][[j]]$k
        
        singlescoreMan[[j]]$dunn <- singlescoreMan[[j]]$dunn + resultMaster[[i]]$hclust$manhattan[[k]][[j]]$dunn
        singlescoreMan[[j]]$sil <- singlescoreMan[[j]]$sil + resultMaster[[i]]$hclust$manhattan[[k]][[j]]$sil
        singlescoreMan[[j]]$gind <- singlescoreMan[[j]]$gind+ resultMaster[[i]]$hclust$manhattan[[k]][[j]]$gind
        singlescoreMan[[j]]$k <- singlescoreMan[[j]]$k + resultMaster[[i]]$hclust$manhattan[[k]][[j]]$k
      }
      singlescore[[j]]$dunn <- singlescore[[j]]$dunn/10
      singlescore[[j]]$sil <- singlescore[[j]]$sil/10
      singlescore[[j]]$gind <- singlescore[[j]]$gind/10
      singlescore[[j]]$k <- singlescore[[j]]$k/10
      
      singlescoreMax[[j]]$dunn <- singlescoreMax[[j]]$dunn/10
      singlescoreMax[[j]]$sil <- singlescoreMax[[j]]$sil/10
      singlescoreMax[[j]]$gind <- singlescoreMax[[j]]$gind/10
      singlescoreMax[[j]]$k <- singlescoreMax[[j]]$k/10
      
      singlescoreMan[[j]]$dunn <- singlescoreMan[[j]]$dunn/10
      singlescoreMan[[j]]$sil <- singlescoreMan[[j]]$sil/10
      singlescoreMan[[j]]$gind <- singlescoreMan[[j]]$gind/10
      singlescoreMan[[j]]$k <- singlescoreMan[[j]]$k/10
    }
    hscores[[i]] <- list(euclidean=singlescore, maximum=singlescoreMax, manhattan=singlescoreMan)
  }
  return(hscores)
  
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

# dla populacji 100 osobnikow, wybieramy co 1/popNr iteracji alg. ewolucyjnego
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

assessGroupingAlgorithm <- function(data, npops, algorithm, func, hmethod = "average", metric = "euclidean", kMin = 2, kMax = 15, tryAllK = FALSE){
  groupingResults = list()
  step = nrow(data)/npops
  
  for(i in 0:(npops-1)){
    start <- i*step+1
    end <- i*step+100
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

populations <- list()
for(funNr in 6:15){
  funcpop <- list()
  for(i in 1:10){
    currentpop <- list()
    rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=100, iters = 1000,
       mutationChance=NA, elitism=NA, monitorFunc=partial(rbgaSavePopulation, name="currentpop"),
                                                        evalFunc = partial(cec2013, i=funNr))
    funcpop[[i]] <- currentpop
  }
  populations[[funNr]] <- funcpop
}

resultMaster <- list()
npops <- 50
hclustMetrics <- list("euclidean", "maximum", "manhattan")
pamMetrics <- list("euclidean", "manhattan")
for(funNr in 6:15){
  resultSingle <- list()

  data_pops <- list()
  #wyciecie npops populacji ze zbioru wygenerowanych danych dla kazdego z 10 przebiegow optymalizacji
  for(n in 1:10){
    data_pops[[n]] <- populationToClusterAnalysis(populationToData(populations[[funNr]][[n]]), npops)
  }
  
  hclust <- list()
  #uruchomienie grupowania dla kazdej z metryk hclust na wycietych populacjach kazdego przebiegu
  for(m in 1:length(hclustMetrics)){
    hclustSingle <- list()
    for(n in 1:10){
      hclustSingle[[n]] <- assessGroupingAlgorithm(data_pops[[n]], npops, "hclust", partial(cec2013, i=funNr), 
                                                   metric = hclustMetrics[[m]])
    }
    hclust[[hclustMetrics[[m]]]] <- hclustSingle
  }
  resultSingle$hclust <- hclust
  
  pamc <- list()
  #uruchomienie grupowania dla kazdej z metryk pam na wycietych populacjach kazdego przebiegu
  for(m in 1:length(pamMetrics)){
    pamcSingle <- list()
    for(n in 1:10){
      pamcSingle[[n]] <- assessGroupingAlgorithm(data_pops[[n]], npops, "pam", partial(cec2013, i=funNr), 
                                                   metric = pamMetrics[[m]])
      pamc[[pamMetrics[[m]]]] <- pamcSingle
    }
    resultSingle$pam <- pamc
  }
  
  kmean <- list()
  #uruchomienie grupowania dla kmeans na wycietych populacjach kazdego przebiegu
  for(n in 1:10){
    kmean[[n]] <- assessGroupingAlgorithm(data_pops[[n]], npops, "kmeans", partial(cec2013, i=funNr))
  }
  resultSingle$kmeans <- kmean
  
  #zrzucenie danych dla danej funkcji cec do glownej listy
  resultMaster[[funNr]] <- resultSingle
}

hscores = list()
for(i in 6:15){
  singlescore = list()
  singlescoreDunn = list()
  singlescoreSil = list()
  singlescoreGind = list()
  singlescoreK = list()
  for(j in seq(1,50)){
    singlescore[[j]] <- list()
    singlescore[[j]]$dunn <- 0
    singlescore[[j]]$sil <- 0
    singlescore[[j]]$gind <- 0
    singlescore[[j]]$k <- 0
    singlescoreDunn <- 0
    singlescoreSil <- 0
    singlescoreGind <- 0
    singlescoreK <- 0 
    for(k in 1:10){
    singlescore[[j]]$dunn <- singlescore[[j]]$dunn + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$dunn
    singlescore[[j]]$sil <- singlescore[[j]]$sil + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$sil
    singlescore[[j]]$gind <- singlescore[[j]]$gind+ resultMaster[[i]]$hclust$euclidean[[k]][[j]]$gind
    singlescore[[j]]$k <- singlescore[[j]]$k + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$k
    singlescoreDunn[[j]] <- singlescoreDunn[[j]] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$dunn
    singlescoreSil[[j]] <- singlescoreSil[[j]] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$sil
    singlescoreGind[[j]] <- singlescoreGind[[j]] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$gind
    singlescoreK[[j]] <- singlescoreK[[j]] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$k
    }
    singlescore[[j]]$dunn <- singlescore[[j]]$dunn/10
    singlescore[[j]]$sil <- singlescore[[j]]$sil/10
    singlescore[[j]]$gind <- singlescore[[j]]$gind/10
    singlescore[[j]]$k <- singlescore[[j]]$k/10
    singlescoreDunn[[j]] <- singlescoreDunn[[j]]/10
    singlescoreSil[[j]] <- singlescoreSil[[j]]/10 
    singlescoreGind[[j]] <- singlescoreGind[[j]]/10
    singlescoreK[[j]] <- singlescoreK[[j]]/10
  }
  hscores[[i]] <- list(singlescoreDunn, singlescoreSil, singlescoreGind, singlescoreK)
}

klist = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("hclust", "pam", "kmeans")))
  for(j in 1:50){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$k
      singlek[j,2] <- singlek[j,2] +resultMaster[[i]]$pam$euclidean[[k]][[j]]$k
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$kmeans[[k]][[j]]$k
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,3]/10
  }
  klist[[i]] <- singlek
}

p1pops <- list()
for(funNr in 6:8){
  p1pops[[funNr]] <- list()
  p1pops[[funNr]]$RBGA <- list()
  p1pops[[funNr]]$GA <- list()
  p1pops[[funNr]]$DE <- list()
  for(i in 1:5){
    currentpop = list()
    rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=100, iters = 1000,
                                      mutationChance=NA, elitism=NA, monitorFunc=partial(rbgaSavePopulation, name="currentpop"),
                                      evalFunc = partial(cec2013, i=funNr))
    p1pops[[funNr]]$RBGA[[i]] <- currentpop
    currentpop = list()
    registerDoSEQ()
    p1pops[[funNr]]$DE[[i]] <- DEoptim(partial(cec2013, i=funNr), rep(-100, 10), rep(100, 10),
                                       DEoptim.control(storepopfrom = 0, trace=FALSE, parallelType=2, itermax=1000))$member$storepop
    currentpop = list()
    ga(type = "real-valued", fitness = partial(cec2013, i=funNr), min = rep(-100, 10), max = rep(100, 10),
       maxiter = 1000, popSize=100, parallel = TRUE, monitor = partial(gaSavePopulation, name="currentpop"))
    p1pops[[funNr]]$GA[[i]] <- currentpop
  }
}

compResults = list()
for(funNr in 6:8){
  compSingle <- list()
  for(alg in c("RBGA", "GA", "DE")){
    for(method in c("hclust", "kmeans"))
    {
      avgSil <- rep(0,10)
      avgDun <- rep(0,10)
      avgGin <- rep(0,10)
      avgK <- rep(0,10)
      for(i in 1:5){
        cur_pops <- populationToClusterAnalysis(populationToData(p1pops[[funNr]][[alg]][[i]]), 10)
        alg_c <- assessGroupingAlgorithm(cur_pops, 10, method, partial(cec2013,funNr))
        for(j in 1:10){
          avgSil[[j]] <- avgSil[[j]] + alg_c[[j]]$sil
          avgDun[[j]] <- avgDun[[j]] + alg_c[[j]]$dunn
          avgK[[j]] <- avgK[[j]] + alg_c[[j]]$k
          avgGin[[j]] <- avgGin[[j]] + alg_c[[j]]$gind
        }
      }
      avgSil <- avgSil/5
      avgDun <- avgDun/5
      avgK <- avgK/5
      avgGin <- avgGin/5
      compSingle[[alg]] <- list()
      compSingle[[alg]][[method]] <- matrix(data = 0, nrow = 10, ncol = 4, dimnames = list(c(), c("sil", "dun", "gin", "k")))
      for(i in 1:10){
        compSingle[[alg]][[method]][i,1] <- avgSil[[i]]
        compSingle[[alg]][[method]][i,2] <- avgDun[[i]]
        compSingle[[alg]][[method]][i,3] <- avgGin[[i]]
        compSingle[[alg]][[method]][i,4] <- avgK[[i]]
      }
    }
  }
  compResults[[funNr]] <- compSingle
}