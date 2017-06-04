plotMetrics<-function(data, numplots, mainlabel, ylabel, maxy){
  x <- seq(1,50,1)
  plot(x, data[1:50,1], col="green", main=mainlabel,
       xlab="nr pokolenia", ylab=ylabel, ylim=c(0,maxy))
  if(numplots==3){
    points(x, data[1:50,2], col="red")
    points(x, data[1:50,3], col="blue")
    legend("topleft", legend=c("euclidean", "manhattan", "maximum"), col=c("green", "red", "blue"),lty=1, cex=1.0)
  }
  else{
    legend("topleft", legend=c("euclidean"), col=c("green"),lty=1, cex=1.0)
  }
}

plotMetrics(data = avghdunn, 3, "Agnes: sredni indeks Dunna dla cec2013(6-15)", "avgDunn", 1.5)
plotMetrics(data = avghsil, 3, "Agnes: sredni wspolczynnik Silhouette dla cec2013(6-15)", "avgSil", 1)
plotMetrics(data = avghgin, 3, "Agnes: sredni odsetek grup o srednich lepszych niz najlepszy dla cec2013(6-15)", "avgGin", 1)
plotMetrics(data = avghk, 3, "Agnes: srednia liczba klastrow dla cec2013(6-15)", "avgSil", 15)

plotMetrics(data = avgpdunn, 3, "PAM: sredni indeks Dunna dla cec2013(6-15)", "avgDunn", 1.5)
plotMetrics(data = avgpsil, 3, "PAM: sredni wspolczynnik Silhouette dla cec2013(6-15)", "avgSil", 1)
plotMetrics(data = avgpgin, 3, "PAM: sredni odsetek grup o srednich lepszych niz najlepszy dla cec2013(6-15)", "avgGin", 1)
plotMetrics(data = avgpk, 3, "PAM: srednia liczba klastrow dla cec2013(6-15)", "avgSil", 15)

plotMetrics(data = avgkdunn, 1, "K-means: sredni indeks Dunna dla cec2013(6-15)", "avgDunn", 1.5)
plotMetrics(data = avgksil, 1, "K-means: sredni wspolczynnik Silhouette dla cec2013(6-15)", "avgSil", 1)
plotMetrics(data = avgkgin, 1, "K-means: sredni odsetek grup o srednich lepszych niz najlepszy dla cec2013(6-15)", "avgGin", 1)
plotMetrics(data = avgkk, 1, "K-means: srednia liczba klastrow dla cec2013(6-15)", "avgSil", 15)