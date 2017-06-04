#Michal Kosikowski, Monika Seniut
#Skrypt do testowania jakosci grupowania populacji tworzonych przez algorytmy ewolucyjne

#Wyniki zebrane w ramach uruchomienia trafiaja do listy zagniezdzonej resultMaster!

#Liczba populacji wycietych rownomiernie sposrod wszystkich wygenerowanych
npops <- 50
#Rozmiar populacji
popsize <- 100
#Liczba iteracji algorytmu RBGA
iterations <- 1000
#Rozpatrywane metryki dla funkcji agnes
agnesMetrics <- list("euclidean", "maximum", "manhattan")
#Rozpatrywane metryki dla funkcji pam
pamMetrics <- list("euclidean", "maximum", "manhattan")
#Parametr method przekazywany do funkcji agnes - okresla metode klastrowania
agnesMethod <- "average"
#Nr pierwszej optymalizowanej funkcji z pakietu cec2013
minFunc <- 6
#Nr ostatniej optymalizowanej funkcji z pakietu cec2013
maxFunc <- 15
#Liczba powtorzen optymalizacji dla kazdej z funkcji
reruns <- 10
#Flaga wlaczajaca grupowanie metoda kmeans
kMeansFlag <- TRUE
#Flaga wlaczajaca grupowanie metoda agnes
agnesFlag <- TRUE
#Flaga wlaczajaca grupowanie metoda pam
pamFlag <- TRUE
#Flaga wlaczajaca generacje grupowania dla wszystkich k z zakresu 2-15
allKFlag <- FALSE

library(cec2013)
library(DEoptim)
library(genalg)
library(GA)
library(cluster)
library(clValid)
library(parallel)
library(doParallel)
library(dbscan)

#Funkcja pomocnicza do tworzenia funkcji czesciowych
partial <- function(f, ...) {
  l <- list(...)
  function(...) {
    do.call(f, c(l, list(...)))
  }
}

#Funkcja pomocnicza do zapisywania populacji stworzonych przez algorytm GA
gaSavePopulation <- function(obj, name){
  env = globalenv()
  env[[name]] = append(env[[name]], list(obj@population))
}

#Funkcja pomocnicza do zapisywania populacji stworzonych przez algorytm RBGA
rbgaSavePopulation <- function(obj, name){
  env = globalenv()
  env[[name]] = append(env[[name]], list(obj$population))
}

populationToData <- function(obj){
  return(do.call(rbind.data.frame, obj))
}

#Zwraca wynik standardowego agnesa z podzialem na k klastrow
#w zakresie minimum:maximum, gdzie k jest wybierane na podstawie
#najwyzszej wartosci sredniej silhouette dla zbioru
#distance - macierz odleglosci grupowanych danych
#type - metoda grupowania sporod dopuszczalnych przez agnes
getBestAgnes <- function(minimum = 2, maximum, distance, type = "average"){
  fit = agnes(distance, method = type)
  bestGrouping = as.data.frame(rep(1, popsize))
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
#dataset - zbior danych jako data frame
#grouping - grupowanie danych jako wektor
#fun - funkcja optymalizowana
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
# wyznaczenie grupowania o najlepszej liczbie dla algorytmu k-srednich
# na podstawie wspolczynnika Silhouette
# k wyznaczane miedzy minimum, a maximum
# distance - macierz odleglosci danych
kMeansClustering <- function(distance, minimum=2, maximum){
  bestGrouping = as.data.frame(rep(1, popsize))
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

# wyznaczenie grupowania o najlepszej liczbie klastrow dla algorytmu PAM
# na podstawie wspolczynnika Silhouette
# k wyznaczane miedzy minimum, a maximum
# distance - macierz odleglosci danych
pamClustering <- function(distance, minimum=2, maximum){
  bestGrouping = as.data.frame(rep(1, popsize))
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

# dla populacji popsize osobnikow, wybieramy co 1/popNr iteracji alg. ewolucyjnego
populationToClusterAnalysis <- function(data, popNr){
  popToClust = matrix(0, nrow = popNr*popsize, ncol = 10)
  popToClust <- data.frame(popToClust)
  coeff = nrow(data)/popNr
  #print("ROWS = ")
  #print(nrow(popToClust))
  #print("Coeff= ")
  #print(coeff)
  for (i in 1:popNr){
    ind1low = (i-1)*popsize+1
    ind1high= ind1low + 99
    ind2low = i*coeff-99 # nie daje (i-1), zeby miec ostatnia populacje
    ind2high = ind2low + 99
    popToClust[ind1low:ind1high,] = data[ind2low:ind2high,]
  }
  # zwraca zbior populacji z poszczegolnych krokow AE, dla ktorej bedzie przeprowadzona analiza grupowania
  return(popToClust)
}

# glowna funkcja zwracajaca wyniki grupowan na danych w zaleznosci od parametrow ustawionych na poczatku pliku
assessGroupingAlgorithm <- function(data, npops, algorithm, func, hmethod = "average", metric = "euclidean", kMin = 2, kMax = 15, tryAllK = FALSE){
  groupingResults = list()
  step = nrow(data)/npops
  
  for(i in 0:(npops-1)){
    start <- i*step+1
    end <- i*step+popsize
    pop.current <- data[c(start:end),]
    if(algorithm != "kmeans" && !is.null(metric)){
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
        if(algorithm == "agnes"){
          kRes$grouping = getBestAgnes(k, k, pop.dist, hmethod)
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
      if(algorithm == "agnes"){
        kRes$grouping = getBestAgnes(kMin, kMax, pop.dist, hmethod)
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

#generacja listy populacji zgodnie z parametrami poczatkowymi skryptu
populations <- list()
for(funNr in minFunc:maxFunc){
  funcpop <- list()
  for(i in 1:reruns){
    currentpop <- list()
    rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=popsize, iters = iterations,
       mutationChance=NA, elitism=NA, monitorFunc=partial(rbgaSavePopulation, name="currentpop"),
                                                        evalFunc = partial(cec2013, i=funNr))
    funcpop[[i]] <- currentpop
  }
  populations[[funNr]] <- funcpop
}

#Lista zagniezdzona zawierajaca komplet wynikow
#struktura:
#resultMaster[[x]]$nazwa_alg_grupowania$metryka[[nr_powtorzenia]][[nr_wycietej_populacji]]$wskaznik
#x - nr funkcji cec
#nazwa_alg_grupowania - agnes, pam lub kmeans
#metryka - dowolna z zawartych w wektorze metryk dla PAM lub agnes,
#ten poziom zagniezdzenia nie jest obecny w przypadku metody kmeans!
#po $kmeans nastepuje od razu[[nr_powtorzenia]]
#nr_powtorzenia - nr z zakresu 1:reruns, okresla nr powtorzenia optymalizacji danej funkcji
#nr_wycietej_populacji - nr z zakresu 1:npops, okresla nr wycietej populacji
#wskaznik - grouping - uzyskane grupowanie danej populacji
#         - dunn - obliczony indeks Dunna danego grupowania
#         - sil - obliczony wspolczynnik Silhouette danego grupowania
#         - gind - obliczony odsetek grup o srednich lepszych niz najlepszy osobnik grupy
#         - k - uzyskana liczba grup
resultMaster <- list()
for(funNr in 6:15){
  resultSingle <- list()

  data_pops <- list()
  #wyciecie npops populacji ze zbioru wygenerowanych danych dla kazdego z 10 przebiegow optymalizacji
  for(n in 1:reruns){
    data_pops[[n]] <- populationToClusterAnalysis(populationToData(populations[[funNr]][[n]]), npops)
  }
  
  if(hClustFlag){
    agnes <- list()
    #uruchomienie grupowania dla kazdej z metryk agnes na wycietych populacjach kazdego przebiegu
    for(m in 1:length(agnesMetrics)){
      agnesSingle <- list()
      for(n in 1:reruns){
        agnesSingle[[n]] <- assessGroupingAlgorithm(data_pops[[n]], npops, "agnes", partial(cec2013, i=funNr),
                                                     metric = agnesMetrics[[m]], hmethod = agnesMethod,
                                                    tryAllK = allKFlag)
      }
      agnes[[agnesMetrics[[m]]]] <- agnesSingle
    }
    resultSingle$agnes <- agnes
  }
  if(pamFlag){
    pamc <- list()
    #uruchomienie grupowania dla kazdej z metryk pam na wycietych populacjach kazdego przebiegu
    for(m in 1:length(pamMetrics)){
      pamcSingle <- list()
      for(n in 1:reruns){
        pamcSingle[[n]] <- assessGroupingAlgorithm(data_pops[[n]], npops, "pam", partial(cec2013, i=funNr), 
                                                     metric = pamMetrics[[m]])
        pamc[[pamMetrics[[m]]]] <- pamcSingle
      }
      resultSingle$pam <- pamc
    }
  }
  
  if(kMeansFlag){
    kmean <- list()
    #uruchomienie grupowania dla kmeans na wycietych populacjach kazdego przebiegu
    for(n in 1:reruns){
      kmean[[n]] <- assessGroupingAlgorithm(data_pops[[n]], npops, "kmeans", partial(cec2013, i=funNr))
    }
    resultSingle$kmeans <- kmean
  }
  #zrzucenie danych dla danej funkcji cec do glownej listy
  resultMaster[[funNr]] <- resultSingle
}