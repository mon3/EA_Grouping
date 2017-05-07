# wczytywanie argumentow

library(cec20005benchmark)
library(cec2013)

population = matrix(nrow = 100, ncol = 10)
# population = apply(population, 1, function(x){runif(100,-100,100)})
 # population = apply(population, 1, function(x){runif(1,-100,100)})
population = replicate(10,runif(100,-100,100))

# goal function evaluation using library 2013 (nr > 5)
cec2013_func_nr = 6
cec2013(cec2013_func_nr, population)

# goal function evaluation using library 2005 (nr > 5)
cec2005benchmark6(population)

# c function return vector
# rep(-10,100) vector of length 100 filled with values: -10


# ga return results for each iteration: iter, population, fitness.
# best, mean, fitnessValue, solution

# TODO: not working with cec fitness functions
# ga(type = "real-valued", fitness =cec2005benchmark6, min=rep(-10,100), max=rep(10,100), nBits = 10)


# genalg RBGA optimization of evaluation function: cec2005benchmark6
library(genalg)
# suggestions == initial population: currently does not work while passing population matrix
rbga.results = rbga(stringMin = rep(-100,10), stringMax = rep(100,10), suggestions=NULL, popSize=100, iters = 100, mutationChance=NA, elitism=NA, monitorFunc=NULL, evalFunc = cec2005benchmark6)
plot(rbga.results)


