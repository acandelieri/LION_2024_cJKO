rm(list=ls()); graphics.off(); cat("\014")

case <- "G2Mx"
N <- 30
run <- 3

source("core.R")

JKO <- readRDS(paste0("JKO_results_N",N,"/JKO_",case,"_FisW2_",run,".RData"))
cJKO <- readRDS(paste0("cJKO_results_N",N,"/cJKO_",case,"_FisW2_",run,".RData"))




Fs <- NULL
Y <- JKO$Y
for( Xk in JKO$Xks ) {
  Fs <- c(Fs,F.sqdW2(Xk))
}

par(mfrow=c(2,2))
par(mar=c(4.1,4.1,1.1,1.1))

plot( 1:length(cJKO$objs), cJKO$objs, type="l", col="blue", lwd=3, ylim=range(Fs,cJKO$objs))
lines( 1:length(Fs), Fs, lwd=3, col="red" )

plot( cumsum(JKO$runTimes), Fs, type="l", col="red", lwd=3,
      ylim=range(Fs,cJKO$objs))
lines( cumsum(cJKO$runTimes), cJKO$objs, lwd=3, col="blue" )


plot( JKO$X0[,1], JKO$X0[,2], pch=19, col="pink", xlim=c(-5,15), ylim=c(-10,10) )
for( i in 1:nrow(JKO$X0) ) {
  lines( c(JKO$X0[i,1],JKO$Xks[[1]][i,1]), c(JKO$X0[i,2],JKO$Xks[[1]][i,2]), col="pink" )
  for( k in 2:length(JKO$Xks) )
    lines( c(JKO$Xks[[k-1]][i,1],JKO$Xks[[k]][i,1]),
           c(JKO$Xks[[k-1]][i,2],JKO$Xks[[k]][i,2]), col="pink")
}
points( JKO$X0[,1], JKO$X0[,2], pch=19, col="red" )
points( JKO$Xks[[length(JKO$Xks)]][,1], JKO$Xks[[length(JKO$Xks)]][,2], pch=4, col="red" )
points( JKO$Y[,1], JKO$Y[,2], pch=1 )


plot( cJKO$X0[,1], cJKO$X0[,2], pch=19, col="deepskyblue", xlim=c(-5,15), ylim=c(-10,10) )
for( i in 1:nrow(cJKO$X0) ) {
  lines( c(cJKO$X0[i,1],cJKO$Xks[[1]][i,1]), c(cJKO$X0[i,2],cJKO$Xks[[1]][i,2]), col="deepskyblue" )
  for( k in 2:length(cJKO$Xks) )
    lines( c(cJKO$Xks[[k-1]][i,1],cJKO$Xks[[k]][i,1]),
           c(cJKO$Xks[[k-1]][i,2],cJKO$Xks[[k]][i,2]), col="deepskyblue")
}
points( cJKO$X0[,1], cJKO$X0[,2], pch=19, col="blue" )
points( cJKO$Xks[[length(cJKO$Xks)]][,1], cJKO$Xks[[length(cJKO$Xks)]][,2], pch=4, col="blue" )
points( cJKO$Y[,1], cJKO$Y[,2], pch=1 )


par(mfrow=c(1,1))
