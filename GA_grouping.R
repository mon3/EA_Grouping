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

RBGA.results = rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=100, iters = 100, 
                    mutationChance=NA, elitism=NA, monitorFunc=partial(rbgaSavePopulation, name="RBGA.pop"),
                    evalFunc = partial(cec2013, i=as.numeric(args[2])))

# RBGA.results = rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=100, iters = 100, 
#                     mutationChance=NA, elitism=NA, monitorFunc=partial(rbgaSavePopulation), name="RBGA.pop", evalFunc = partial(cec2013, i=7)) # i to make a parameter
# partial(cec2013, i=7) zwraca jako wynik jedną liczbę dla każdego wektora (osobnika) z populacji

# plot(RBGA.results)



# run before parallel DE to avoid errors
registerDoSEQ()
DE.results = DEoptim(partial(cec2013, i=as.numeric(args[2])), rep(-100, 10), rep(100, 10), DEoptim.control(storepopfrom = 0, trace=FALSE, parallelType=2))
# plot(DEres, plot.type = "bestvalit")
DE.pop = DE.results$member$storepop # wycinamy member storepop
GA.results = ga(type = "real-valued", fitness = partial(cec2013, i=as.numeric(args[2])), min = rep(-100, 10), max = rep(100, 10),
           maxiter = 500, popSize=100, parallel = TRUE, monitor = partial(gaSavePopulation, name="GA.pop"))
plot(GA.results)


# distance matrix
# d <- dist(mydata, method = "euclidean") # distance matrix
