hdunnscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$agnes$euclidean[[k]][[j]]$dunn
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$agnes$manhattan[[k]][[j]]$dunn
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$agnes$maximum[[k]][[j]]$dunn
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,3]/10
  }
  hdunnscores[[i]] <- singlek
}

avghdunn <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
for(i in 6:15){
  avghdunn <- avghdunn + hdunnscores[[i]]
}
avghdunn <- avghdunn/10


hsilscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$agnes$euclidean[[k]][[j]]$sil
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$agnes$manhattan[[k]][[j]]$sil
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$agnes$maximum[[k]][[j]]$sil
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,3]/10
  }
  hsilscores[[i]] <- singlek
}

avghsil <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
for(i in 6:15){
  avghsil <- avghsil + hsilscores[[i]]
}
avghsil <- avghsil/10

hginscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$agnes$euclidean[[k]][[j]]$gind
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$agnes$manhattan[[k]][[j]]$gind
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$agnes$maximum[[k]][[j]]$gind
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,3]/10
  }
  hginscores[[i]] <- singlek
}

avghgin <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
for(i in 6:15){
  avghgin <- avghgin + hginscores[[i]]
}
avghgin <- avghgin/10

hkscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$agnes$euclidean[[k]][[j]]$k
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$agnes$manhattan[[k]][[j]]$k
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$agnes$maximum[[k]][[j]]$k
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,3]/10
  }
  hkscores[[i]] <- singlek
}

avghk <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
for(i in 6:15){
  avghk <- avghk + hkscores[[i]]
}
avghk <- avghk/10


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
    singlek[j,3] <- singlek[j,3]/10
  }
  pdunnscores[[i]] <- singlek
}

avgpdunn <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
for(i in 6:15){
  avgpdunn <- avgpdunn + pdunnscores[[i]]
}
avgpdunn <- avgpdunn/10

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
    singlek[j,3] <- singlek[j,3]/10
  }
  psilscores[[i]] <- singlek
}

avgpsil <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
for(i in 6:15){
  avgpsil <- avgpsil + psilscores[[i]]
}
avgpsil <- avgpsil/10

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
    singlek[j,3] <- singlek[j,3]/10
  }
  pginscores[[i]] <- singlek
}

avgpgin <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
for(i in 6:15){
  avgpgin <- avgpgin + pginscores[[i]]
}
avgpgin <- avgpgin/10

pkscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$pam$euclidean[[k]][[j]]$k
      singlek[j,2] <- singlek[j,2] + resultMaster[[i]]$pam$manhattan[[k]][[j]]$k
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$pam$maximum[[k]][[j]]$k
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,3]/10
  }
  pkscores[[i]] <- singlek
}

avgpk <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("euclidean", "manhattan", "maximum")))
for(i in 6:15){
  avgpk <- avgpk + pkscores[[i]]
}
avgpk <- avgpk/10

kdunnscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$kmeans[[k]][[j]]$dunn
    }
    singlek[j,1] <- singlek[j,1]/10
  }
  kdunnscores[[i]] <- singlek
}

avgkdunn <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
for(i in 6:15){
  avgkdunn <- avgkdunn + kdunnscores[[i]]
}
avgkdunn <- avgkdunn/10

ksilscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$kmeans[[k]][[j]]$sil
    }
    singlek[j,1] <- singlek[j,1]/10
  }
  ksilscores[[i]] <- singlek
}

avgksil <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
for(i in 6:15){
  avgksil <- avgksil + ksilscores[[i]]
}
avgksil <- avgksil/10

kginscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$kmeans[[k]][[j]]$gind
    }
    singlek[j,1] <- singlek[j,1]/10
  }
  kginscores[[i]] <- singlek
}

avgkgin <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
for(i in 6:15){
  avgkgin <- avgkgin + kginscores[[i]]
}
avgkgin <- avgkgin/10

kkscores = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
  for(j in seq(1,50)){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$kmeans[[k]][[j]]$k
    }
    singlek[j,1] <- singlek[j,1]/10
  }
  kkscores[[i]] <- singlek
}

avgkk <- matrix(data = 0, nrow = 50, ncol = 1, dimnames = list(c(), c("euclidean")))
for(i in 6:15){
  avgkk <- avgkk + kkscores[[i]]
}
avgkk <- avgkk/10

klist = list()
for(i in 6:15){
  singlek <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("agnes", "pam", "kmeans")))
  for(j in 1:50){
    for(k in 1:10){
      singlek[j,1] <- singlek[j,1] + resultMaster[[i]]$agnes$euclidean[[k]][[j]]$k
      singlek[j,2] <- singlek[j,2] +resultMaster[[i]]$pam$euclidean[[k]][[j]]$k
      singlek[j,3] <- singlek[j,3] + resultMaster[[i]]$kmeans[[k]][[j]]$k
    }
    singlek[j,1] <- singlek[j,1]/10
    singlek[j,2] <- singlek[j,2]/10
    singlek[j,3] <- singlek[j,3]/10
  }
  klist[[i]] <- singlek
}

avgbestdunn <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("kmeans", "agnes euclidean", "pam maximum")))
avgbestdunn[,1] <- avgkdunn[,1]
avgbestdunn[,2] <- avghdunn[,1]
avgbestdunn[,3] <- avgpdunn[,3]

avgbestsil <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("kmeans", "agnes euclidean", "pam maximum")))
avgbestsil[,1] <- avgksil[,1]
avgbestsil[,2] <- avghsil[,1]
avgbestsil[,3] <- avgpsil[,3]

avgbestgin <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("kmeans", "agnes euclidean", "pam maximum")))
avgbestgin[,1] <- avgkgin[,1]
avgbestgin[,2] <- avghgin[,1]
avgbestgin[,3] <- avgpgin[,3]

avgbestk <- matrix(data = 0, nrow = 50, ncol = 3, dimnames = list(c(), c("kmeans", "agnes euclidean", "pam maximum")))
avgbestk[,1] <- avgkk[,1]
avgbestk[,2] <- avghk[,1]
avgbestk[,3] <- avgpk[,3]