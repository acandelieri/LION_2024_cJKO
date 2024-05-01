rm(list=ls()); graphics.off(); cat("\014")

source("core.R")

d <- 2
N <- 15

K.jko <- 1000 
patience.jko <- 10
h <- 0.01

K.cjko <- 5000 
patience.cjko <- 10
eps <- 0.01




Ff <- F.sqdW2_fromGaussian
set.seed(42)
GY <- rmvnorm( n=N, mean=c(3,3), sigma=0.1*diag(x=1,nrow=d) )


X0 <- rmvnorm( n=N, mean=c(-3,-3), sigma=0.2*diag(1,nrow=d) )

cat("> Starting JKO\n")
jko.res <- JKO( X=X0, Ff=Ff, h=h, K=K.jko, patience=patience.jko )
cat("\n> Solved in:",jko.res$elapsed,"[secs]\n")

cat("\n> Starting constrained-JKO\n")
cjko.res <- cJKO( X=X0, Ff=Ff, eps=eps, K=K.cjko, patience=patience.cjko )
cat("\n> Solved in:",cjko.res$elapsed,"[secs]\n")





plot( jko.res$objs, type="l", ylim=range(jko.res$objs),
      lwd=3, xlab="iterations", ylab="0.5 W_2^2 + F", col="blue" )
cjko.cmp <- cjko.res$W2s/(2*h) + cjko.res$objs
lines( cjko.cmp, lwd=3, col="green3" )


plot( cjko.res$objs, type="l", lwd=3, xlab="cJKO iterations", ylab="F(P)", col="green3" )

stop("FINE")

# for( k in 1:K ) {
#   if( k==1 ) {
#     plot( jko.res$Xks[[k]][,1], jko.res$Xks[[k]][,2], pch=19, col="blue",
#           xlim=c(-5,5), ylim=c(-5,5))
#   } else {
#     points( jko.res$Xks[[k-1]][,1], jko.res$Xks[[k-1]][,2], pch=19, col="grey" )
#     points( jko.res$Xks[[k]][,1], jko.res$Xks[[k]][,2], pch=19, col="blue" )
#     for( i in 1:N )
#       lines( c(jko.res$Xks[[k-1]][i,1],jko.res$Xks[[k]][i,1]),
#                c(jko.res$Xks[[k-1]][i,2],jko.res$Xks[[k]][i,2]), col="deepskyblue" )
#     Sys.sleep(1)
#   }
# }
# 
# stop("NEXT?")
# 
# for( k in 1:K ) {
#   if( k==1 ) {
#     plot( cjko.res$Xks[[k]][,1], cjko.res$Xks[[k]][,2], pch=19, col="green3",
#           xlim=c(-5,5), ylim=c(-5,5))
#   } else {
#     points( cjko.res$Xks[[k-1]][,1], cjko.res$Xks[[k-1]][,2], pch=19, col="grey" )
#     points( cjko.res$Xks[[k]][,1], cjko.res$Xks[[k]][,2], pch=19, col="green3" )
#     Sys.sleep(1)
#   }
# }
# 
# 
# stop("NEXT?")

plot( X0[,1], X0[,2], pch=19, col="blue",
      xlim=c(-5,5), ylim=c(-5,5))
for( k in 2:length(jko.res$Xks) ) {
  for( i in 1:N )
    lines( x = c(jko.res$Xks[[k-1]][i,1],jko.res$Xks[[k]][i,1]),
          y = c(jko.res$Xks[[k-1]][i,2],jko.res$Xks[[k]][i,2]), col="deepskyblue" )
}
points( jko.res$Xks[[1]][,1], jko.res$Xks[[1]][,2], pch=19, col="blue" )
points( jko.res$Xks[[K.jko]][,1], jko.res$Xks[[K.jko]][,2], pch=8, col="blue" )

points(GY[,1],GY[,2],pch=1.2,lwd=2,col="red")


stop("NEXT?")

plot( X0[,1], X0[,2], pch=19, col="green3",
      xlim=c(-5,5), ylim=c(-5,5))
for( k in 2:length(cjko.res$Xks) ) {
  for( i in 1:N )
    lines( x = c(cjko.res$Xks[[k-1]][i,1],cjko.res$Xks[[k]][i,1]),
           y = c(cjko.res$Xks[[k-1]][i,2],cjko.res$Xks[[k]][i,2]), col="lightgreen" )
}
points( cjko.res$Xks[[1]][,1], cjko.res$Xks[[1]][,2], pch=19, col="green3" )
points( cjko.res$Xks[[K.cjko]][,1], cjko.res$Xks[[K.cjko]][,2], pch=8, col="green3" )

points(GY[,1],GY[,2],pch=1.2,lwd=2,col="red")

