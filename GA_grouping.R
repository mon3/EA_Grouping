# wczytywanie argumentow

library(cec2005benchmark)
library(cec2013)
library(DEoptim)
library(genalg)
library(GA)

library(parallel)
library(doParallel)


# arguments to control included
# arg1 = number of cec2013 function for testing purposes
# arg2 = number of cec2013 used in genetic algorithms using different packages
args <- commandArgs(trailingOnly = TRUE)
print(args[1])
print(args[2])


population = matrix(nrow = 100, ncol = 10)
# population = apply(population, 1, function(x){runif(100,-100,100)})
 # population = apply(population, 1, function(x){runif(1,-100,100)})
population = replicate(10,runif(100,-100,100))

# goal function evaluation using library 2013 (nr > 5)
cec2013_func_nr = as.numeric(args[1]) # 6 - first arguemnt or any number
cec2013(cec2013_func_nr, population) # returns calculated goal function for the population
# 1 number for individual

# generates a partial function - used to pass an arbitrary cec2013 function to DEoptim
partial <- function(f, ...) {
  l <- list(...)
  function(...) {
    do.call(f, c(l, list(...)))
  }
}

#saves intermediate GA populations to an global environment variable named name
gaSavePopulation <- function(obj, name){
  env = globalenv()
  env[[name]] = append(env[[name]], list(obj@population))
}

#saves intermediate RBGA populations to an global environment variable named name
rbgaSavePopulation <- function(obj, name){
  env = globalenv()
  env[[name]] = append(env[[name]], list(obj$population))
}

# revrites set of populations from genetic algorithms to data frame
populationToData <- function(obj){
	return(do.call(rbind.data.frame, obj))
}

# wss <- (nrow(RBGA.lastpop)-1)*sum(apply(RBGA.lastpop,2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(RBGA.lastpop, 
#    centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#    ylab="Within groups sum of squares")


# determination of optimal number of clusters
determineNumberOfClusters <- function(obj){
	wss <- (nrow(obj)-1)*sum(apply(obj,2,var))
	for (i in 2:15) wss[i] <- sum(kmeans(obj, 
   		centers=i)$withinss)
	plot(1:15, wss, type="b", xlab="Number of Clusters",
   		ylab="Within groups sum of squares")
}


# K-means clusters analysis
kMeansClustering <- function(obj, numberOfClusters){
	determineNumberOfClusters(obj)
	# k-means cluster analysis
	fit <- kmeans(obj, numberOfClusters)
	# get cluster means 
	aggregate(obj, by=list(fit$cluster), FUN=mean)
	# append cluster assignment
	obj <- data.frame(obj, fit$cluster)
}


# np. maxiter = 1000, dalej potrzebujemy większej liczby grup niż w przypadku RBGA






# goal function evaluation using library 2005 (nr > 5)
cec2005benchmark6(population)

# results that are easier to interpret
cec2005benchmark9(population)


# c function return vector
# rep(-10,100) vector of length 100 filled with values: -10


# ga return results for each iteration: iter, population, fitness.
# best, mean, fitnessValue, solution

#empty lists for populations
GA.pop = list()
RBGA.pop = list()
# TODO: not working with cec fitness functions
# ga(type = "real-valued", fitness =cec2005benchmark6, min=rep(-10,100), max=rep(10,100), nBits = 10)

# genalg RBGA optimization of evaluation function: cec2005benchmark6
# suggestions == initial population: currently does not work while passing population matrix


# zmieniałam maxiter, żeby sprawdzić, jak się nasyca funkcja celu i jednocześnie, ile grup powstanie
# maxiter = 100, #grup = 14
# maxiter = 500, #grup = 3 (tak naprawdę 1)
# maxiter = 300, #grup = 6/7 

RBGA.results = rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=100, iters = 300, 
                    mutationChance=NA, elitism=NA, monitorFunc=partial(rbgaSavePopulation, name="RBGA.pop"),
                    evalFunc = partial(cec2013, i=9))

# RBGA.results = rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=100, iters = 100, 
#                     mutationChance=NA, elitism=NA, monitorFunc=partial(rbgaSavePopulation), name="RBGA.pop", evalFunc = partial(cec2013, i=7)) # i to make a parameter
# partial(cec2013, i=7) zwraca jako wynik jedną liczbę dla każdego wektora (osobnika) z populacji

# plot(RBGA.results)



# run before parallel DE to avoid errors
registerDoSEQ()
DE.results = DEoptim(partial(cec2013, i=as.numeric(args[2])), rep(-100, 10), rep(100, 10), DEoptim.control(storepopfrom = 0, trace=FALSE, parallelType=2))
# plot(DEres, plot.type = "bestvalit")
DE.pop = DE.results$member$storepop # wycinamy member storepop

GA.results = ga(type = "real-valued", fitness = partial(cec2013, i=7), min = rep(-100, 10), max = rep(100, 10),
           maxiter = 1000, popSize=100, parallel = TRUE, monitor = partial(gaSavePopulation, name="GA.pop"))
plot(GA.results)


# distance matrix
# d <- dist(mydata, method = "euclidean") # distance matrix


RBGA.data = populationToData(RBGA.pop)
RBGA.lastpop = tail(RBGA.data, 100)

GA.data = populationToData(GA.pop)
GA.lastpop = tail(GA.data, 100)

# determine number of clusters for RBGA last population
determineNumberOfClusters(RBGA.lastpop)

# determine number of clusters for GA last population
# z wykresów wynika, że dla GA nawet przy większej liczbie iteracji, 
# np. maxiter = 1000, dalej potrzebujemy większej liczby grup niż w przypadku RBGA
determineNumberOfClusters(GA.lastpop)

# K-Means Cluster Analysis
 fit <- kmeans(RBGA.lastpop, 14) # 6 cluster solution
# get cluster means 
aggregate(RBGA.lastpop, by=list(fit$cluster), FUN =mean)
# append cluster assignment
RBGA.lastpop <- data.frame(RBGA.lastpop, fit$cluster)

# sprawdzenie wynikow grupowania: plot(fit$cluster)


# creation of clusters from k-means for GA
 fit <- kmeans(GA.lastpop, 8) # 8 cluster solution
# get cluster means 
aggregate(GA.lastpop, by=list(fit$cluster), FUN =mean)
# append cluster assignment
GA.lastpop <- data.frame(GA.lastpop, fit$cluster)



# determineNumberOfClusters(GA.lastpop)

# # K-Means Cluster Analysis
#  fit <- kmeans(RBGA.lastpop, 6) # 6 cluster solution
# # get cluster means 
# aggregate(RBGA.lastpop, by=list(fit$cluster), FUN =mean)
# # append cluster assignment
# RBGA.lastpop <- data.frame(RBGA.lastpop, fit$cluster)

# kMeansClustering(RBGA.lastpop, 6)



# ANOTHER POSSIBILITY for k-means
# nstart = number of initial points
#RBGA.Cluster <- kmeans(RBGA.lastpop, 7, nstart =20)


# porównanie funkcji celu dla środków każdej z grup z funkcją celu 
# dla najlepszego osobnika z każdej z grup

# grupa o indeksie 1
RBGA.lastpop.sub3 <-subset(RBGA.lastpop, fit.cluster==3) # zrobić pętlę for
cec2013(9, fit$centers)
group1 = data.matrix(RBGA.lastpop.sub1[,1:10])
group3 = data.matrix(RBGA.lastpop.sub3[,1:10])

min(cec2013(9,group1)) # daje wartośc minimalną: działa dla dodatnich i ujemnych wartości

goalCluster=(cec2013(9,fit$centers))
minPopGoal =min(cec2013(9,group1))
minPopGoal =min(cec2013(9,group3))


# dla ujemnych liczb: porównanie funkcji celu
res = goalCluster>=minPopGoal # muszą być wszystkie TRUE
all(res)==TRUE # sprawdzamy, czy wszystkie są true: rozpatrywany przypadek dla ujemnych liczb!


# podział RBGA.lastpop na podgrupy; RBGA.lastpop.subsets = lista
RBGA.lastpop.subsets <- split(RBGA.lastpop, as.factor(RBGA.lastpop$fit.cluster))


	
cec2013(9,data.matrix(as.data.frame(RBGA.lastpop.subsets[1])[,1:10]))

# wypisuje funckcje celu dla osobników w poszczególnych grupach ->
# należy je porównać z wynikami dla średnich z clusters goalCluster=(cec2013(9,fit$centers)) w grupach
for (i in 1:length(unique(RBGA.lastpop$fit.cluster)))
{
	print(cec2013(9, data.matrix(as.data.frame(RBGA.lastpop.subsets[i])[,1:10])))
}




