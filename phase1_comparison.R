library(cec2013)
library(DEoptim)
library(genalg)
library(GA)
library(cluster)
library(clValid)
library(parallel)
library(doParallel)

#Michal Kosikowski, Monika Seniut
#Pomocniczy skrypt sluzacy do generacji danych do pierwszego etapu badan

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
    for(method in c("kmeans"))
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