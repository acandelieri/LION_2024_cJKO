rm(list=ls()); graphics.off(); cat("\014")

source("core.R")

case <- "G2G"
N <- 15
run <- 5
obj.fun <- "FisW2"
ng <- 30

jko <- readRDS(paste0("JKO_results_N",N,"/JKO_",case,"_",obj.fun,"_",run,".RData"))
cjko <- readRDS(paste0("cJKO_results_N",N,"/cJKO_",case,"_",obj.fun,"_",run,".RData"))

stopifnot( all(jko$Y == cjko$Y) && all(jko$X0 == cjko$X0) )


cat("[ *** Results for",case,"run",run,"*** ]\n\n")

cat("JKO's h =",jko$h,"\n")
cat("cJKO's epsilon =",cjko$eps,"\n\n")


# *******************************************************************
# Monge Gap
# *******************************************************************

jko.cost <- sqrt(sum((jko$X0 - jko$Xks[[length(jko$Xks)]])^2)/nrow(jko$X0))
ot.jko <- transport(pp(jko$X0),pp(jko$Xks[[length(jko$Xks)]]),p=2,method="primaldual")
ot.jko.cost <- wasserstein(pp(jko$X0),pp(jko$Xks[[length(jko$Xks)]]),p=2,tplan=ot.jko )
MG.jko <- round(jko.cost - ot.jko.cost,4)
cat("Monge Gap JKO:",MG.jko,"\n")

cjko.cost <- sqrt(sum((cjko$X0 - cjko$Xks[[length(cjko$Xks)]])^2)/nrow(cjko$X0))
ot.cjko <- transport(pp(cjko$X0),pp(cjko$Xks[[length(cjko$Xks)]]),p=2,method="primaldual")
ot.cjko.cost <- wasserstein(pp(cjko$X0),pp(cjko$Xks[[length(cjko$Xks)]]),p=2,tplan=ot.cjko )
MG.cjko <- round(cjko.cost - ot.cjko.cost,4)
cat("Monge Gap cJKO:",MG.cjko,"\n\n")


# *******************************************************************
# F(P_K)
# *******************************************************************

stopifnot( all(jko$Y==cjko$Y) )
Y <- jko$Y

if( obj.fun=="FisW2" ) {
  F.jko <- F.sqdW2( jko$Xks[[length(jko$Xks)]] )
  F.cjko <- F.sqdW2( cjko$Xks[[length(cjko$Xks)]] )
} else { 
  target.mean <- jko$tgtMeans
  target.sigma <- jko$tgtSigmas
  F.jko <- F.symKL( jko$Xks[[length(jko$Xks)]] )
  F.cjko <- F.symKL( cjko$Xks[[length(cjko$Xks)]] )
}
  
cat("F(P_K) for JKO:",round(F.jko,4),"\n")
cat("F(P_K) for cJKO:",round(F.cjko,4),"\n\n")


# *******************************************************************
# Mismatch (i.e. FID)
# *******************************************************************

ot.jko <- transport(pp(jko$Xks[[length(jko$Xks)]]),pp(jko$Y),p=2,method="primaldual")
mis.jko <- wasserstein(pp(jko$Xks[[length(jko$Xks)]]),pp(jko$Y),p=2,tplan=ot.jko)
cat("Mismatch JKO:",round(mis.jko,4),"\n")

ot.cjko <- transport(pp(cjko$Xks[[length(cjko$Xks)]]),pp(cjko$Y),p=2,method="primaldual")
mis.cjko <- wasserstein(pp(cjko$Xks[[length(cjko$Xks)]]),pp(cjko$Y),p=2,tplan=ot.cjko)
cat("Mismatch cJKO:",round(mis.cjko,4),"\n\n")



# *******************************************************************
# Transport cost (average per iter)
# *******************************************************************
 
jko.Tcost <- 0
for( i in 2:length(jko$Xks) )
   jko.Tcost <- jko.Tcost + sum( sqrt(apply( (jko$Xks[[i-1]] - jko$Xks[[i]])^2, 1 , sum )) )

cjko.Tcost <- 0
for( i in 2:length(cjko$Xks) )
  cjko.Tcost <- cjko.Tcost + sum( sqrt(apply( (cjko$Xks[[i-1]] - cjko$Xks[[i]])^2, 1 , sum )) )

cat("JKO's T cost (avg per iter):",round(jko.Tcost/length(jko$Xks),4),"\n")
cat("cJKO's T cost (avg per iter):",round(cjko.Tcost/length(cjko$Xks),4),"\n\n")


# *******************************************************************
# Number of iterations
# *******************************************************************

cat("JKO iterations:",length(jko$Xks),"\n")
cat("cJKO iterations:", length(cjko$Xks),"\n\n")


# *******************************************************************
# Runtime
# *******************************************************************

cat("Runtime JKO [secs]:", sum(jko$runTimes),"\n")
cat("Runtime cJKO [secs]:", sum(cjko$runTimes),"\n\n")



# *******************************************************************
# Some charts
# *******************************************************************

OT <- transport(pp(jko$X0),pp(jko$Y),p=2,method="primaldual")

par(mar=c(4.1,4.6,3.1,1.1) )
plot( jko$Y[,1], jko$Y[,2], pch=1, cex=1.3, xlim=c(-5,15), ylim=c(-5,10),
      xlab=expression(x[1]), ylab=expression(x[2]), main="JKO",
      cex.main=2, cex.lab=2, cex.axis=2 )
for( i in 1:nrow(jko$X0) ) {
  lines( c(jko$X0[OT$from[i],1], jko$Y[OT$to[i],1]), c(jko$X0[OT$from[i],2], jko$Y[OT$to[i],2]),
         col=adjustcolor("black",alpha.f=0.3), lty=2, lwd=2 )
  lines( c(jko$X0[i,1], jko$Xks[[length(jko$Xks)]][i,1]),
         c(jko$X0[i,2], jko$Xks[[length(jko$Xks)]][i,2]),
         lwd=2, col=adjustcolor("red",alpha.f=0.5) )
}
points( jko$Y[,1], jko$Y[,2], pch=1, cex=1.3 )
points( jko$X0[,1], jko$X0[,2], pch=19, col="red" )
points( jko$Xks[[length(jko$Xks)]][,1], jko$Xks[[length(jko$Xks)]][,2],
        pch=20, col="red", lwd=2 )
legend("bottomright", legend=paste0("h=",jko$h), cex=1.5)



plot( cjko$Y[,1], cjko$Y[,2], pch=1, cex=1.3, xlim=c(-5,15), ylim=c(-5,10),
      xlab=expression(x[1]), ylab=expression(x[2]), main="cJKO",
      cex.main=2, cex.lab=2, cex.axis=2 )
for( i in 1:nrow(cjko$X0) ) {
  lines( c(cjko$X0[OT$from[i],1], cjko$Y[OT$to[i],1]),
         c(cjko$X0[OT$from[i],2], cjko$Y[OT$to[i],2]),
         col=adjustcolor("black",alpha.f=0.3), lty=2, lwd=2 )
  lines( c(cjko$X0[i,1], cjko$Xks[[length(cjko$Xks)]][i,1]),
         c(cjko$X0[i,2], cjko$Xks[[length(cjko$Xks)]][i,2]),
         lwd=2, col=adjustcolor("blue",alpha.f=0.5) )
}
points( cjko$Y[,1], cjko$Y[,2], pch=1, cex=1.3 )
points( cjko$X0[,1], cjko$X0[,2], pch=19, col="blue" )
points( cjko$Xks[[length(cjko$Xks)]][,1], cjko$Xks[[length(cjko$Xks)]][,2],
        pch=20, col="blue", lwd=2 )
legend("bottomright", legend=paste0("eps=",cjko$eps), cex=1.5)


