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

cutAverageResPam <- function(resultMaster){
  hscores = list()
  
  for(i in 6:15){
    singlescore = list()
    singlescoreMan = list()
    for(j in seq(1,50)){
      singlescore[[j]] <- list()
      singlescore[[j]]$dunn <- 0
      singlescore[[j]]$sil <- 0
      singlescore[[j]]$gind <- 0
      singlescore[[j]]$k <- 0
  
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
        
        singlescoreMan[[j]]$dunn <- singlescoreMan[[j]]$dunn + resultMaster[[i]]$hclust$manhattan[[k]][[j]]$dunn
        singlescoreMan[[j]]$sil <- singlescoreMan[[j]]$sil + resultMaster[[i]]$hclust$manhattan[[k]][[j]]$sil
        singlescoreMan[[j]]$gind <- singlescoreMan[[j]]$gind+ resultMaster[[i]]$hclust$manhattan[[k]][[j]]$gind
        singlescoreMan[[j]]$k <- singlescoreMan[[j]]$k + resultMaster[[i]]$hclust$manhattan[[k]][[j]]$k
      }
      singlescore[[j]]$dunn <- singlescore[[j]]$dunn/10
      singlescore[[j]]$sil <- singlescore[[j]]$sil/10
      singlescore[[j]]$gind <- singlescore[[j]]$gind/10
      singlescore[[j]]$k <- singlescore[[j]]$k/10
      
      singlescoreMan[[j]]$dunn <- singlescoreMan[[j]]$dunn/10
      singlescoreMan[[j]]$sil <- singlescoreMan[[j]]$sil/10
      singlescoreMan[[j]]$gind <- singlescoreMan[[j]]$gind/10
      singlescoreMan[[j]]$k <- singlescoreMan[[j]]$k/10
    }
    hscores[[i]] <- list(euclidean=singlescore, manhattan=singlescoreMan)
  }
  return(hscores)
  
}


cutAverageResKmeans <- function(resultMaster){
  hscores = list()
  
  for(i in 6:15){
    singlescore = list()
    singlescoreMan = list()
    for(j in seq(1,50)){
      singlescore[[j]] <- list()
      singlescore[[j]]$dunn <- 0
      singlescore[[j]]$sil <- 0
      singlescore[[j]]$gind <- 0
      singlescore[[j]]$k <- 0
      
      for(k in 1:10){
        singlescore[[j]]$dunn <- singlescore[[j]]$dunn + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$dunn
        singlescore[[j]]$sil <- singlescore[[j]]$sil + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$sil
        singlescore[[j]]$gind <- singlescore[[j]]$gind+ resultMaster[[i]]$hclust$euclidean[[k]][[j]]$gind
        singlescore[[j]]$k <- singlescore[[j]]$k + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$k

      }
      singlescore[[j]]$dunn <- singlescore[[j]]$dunn/10
      singlescore[[j]]$sil <- singlescore[[j]]$sil/10
      singlescore[[j]]$gind <- singlescore[[j]]$gind/10
      singlescore[[j]]$k <- singlescore[[j]]$k/10
    }
    hscores[[i]] <- list(euclidean=singlescore)
  }
  return(hscores)
  
}
