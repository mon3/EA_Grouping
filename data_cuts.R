hdunnscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$dunn
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$hclust$manhattan[[k]][[j]]$dunn
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$hclust$maximum[[k]][[j]]$dunn
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,3]/10
  }
  hdunnscores[[i]] <- singlek
}

hsilscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$sil
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$hclust$manhattan[[k]][[j]]$sil
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$hclust$maximum[[k]][[j]]$sil
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,3]/10
  }
  hsilscores[[i]] <- singlek
}

hginscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$gind
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$hclust$manhattan[[k]][[j]]$gind
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$hclust$maximum[[k]][[j]]$gind
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,3]/10
  }
  hginscores[[i]] <- singlek
}

hkscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$k
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$hclust$manhattan[[k]][[j]]$k
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$hclust$maximum[[k]][[j]]$k
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,3]/10
  }
  hkscores[[i]] <- singlek
}

pdunnscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$pam$euclidean[[k]][[j]]$dunn
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$pam$manhattan[[k]][[j]]$dunn
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$pam$maximum[[k]][[j]]$dunn
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,2]/10
  }
  pdunnscores[[i]] <- singlek
}

psilscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$pam$euclidean[[k]][[j]]$sil
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$pam$manhattan[[k]][[j]]$sil
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$pam$maximum[[k]][[j]]$sil
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,2]/10
  }
  psilscores[[i]] <- singlek
}

pginscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$pam$euclidean[[k]][[j]]$gind
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$pam$manhattan[[k]][[j]]$gind
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$pam$maximum[[k]][[j]]$gind
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,2]/10
  }
  pginscores[[i]] <- singlek
}

pkscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 2, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$pam$euclidean[[k]][[j]]$k
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$pam$manhattan[[k]][[j]]$k
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$pam$maximum[[k]][[j]]$k
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,2]/10
  }
  pkscores[[i]] <- singlek
}

kdunnscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$dunn
    }
    singlek[j,1] <- singlek[j,1]/10
  }
  kkscores[[i]] <- singlek
}

ksilscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$sil
    }
    singlek[j,1] <- singlek[j,1]/10
  }
  kkscores[[i]] <- singlek
}

kginscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$gind
    }
    singlek[j,1] <- singlek[j,1]/10
  }
  kkscores[[i]] <- singlek
}

kkscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$hclust$euclidean[[k]][[j]]$k
    }
    singlek[j,1] <- singlek[j,1]/10
  }
  kkscores[[i]] <- singlek
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