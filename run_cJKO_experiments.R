rm(list=ls()); graphics.off(); cat("\014")

source("core.R")

#*********************************************************************************
# Experiment set up
#*********************************************************************************

d <- 2  # dimensionality of the support
N <- 30 # size of the point cloud

K.cjko <- 5000  # max cJKO iterations
patience.cjko <- 10 # number of consecutive iterations without any improvement

# cJKO's parameter: in experiments eps in {0.5, 0.25, 0.1, 0.05, 0.01} 
eps <- 0.01

# case study
case <- "Mx2G"  # in the experiments one among:
                # 'G2G'  (Gaussian multi-variate distribution to another ),
                # 'G2Mx' (Gaussian multi-variate distribution to a Gaussian Mixture),
                # 'Mx2G' (Gaussian Mixture to a Gaussian multi-variate distribution)

# objective function F(P)
obj.fun <- "FisSymKL" # in the experiments, two possibilities:
                      # 'FisW2', F(P) is the 2-Wasserstein distance from target distribution
                      # 'FisSymKL' F(P) is the symmetrised Kullback-Liebler from the target distribution
ng <- 30 # number of points per dimension for the computation of the 'symKL'





#*********************************************************************************
# MAIN
#*********************************************************************************

if( case == "G2G" ) {
  stopifnot( obj.fun %in% c("FisW2","FisSymKL") )
  Ff <- ifelse(obj.fun=="FisW2",F.sqdW2,F.symKL)
  set.seed(42)
  target.mean <- c(10,5)
  target.sigma <- 2*matrix( c(1,-0.75,-0.75,1), nrow=d )
  Y <- rmvnorm( n=N, mean=target.mean, sigma=target.sigma )
  source.mean <- c(0,0)
  source.sigma <- diag(x=1,nrow=d)
  X0 <- rmvnorm( n=N, mean=source.mean, sigma=source.sigma )
} else {
  if( case == "G2Mx" ) {
    stopifnot(N%%3==0)
    stopifnot( obj.fun %in% c("FisW2","FisSymKL") )
    Ff <- ifelse(obj.fun=="FisW2",F.sqdW2,F.symKL)
    set.seed(42)
    target.mean <- list( c(10,7), c(5,5), c(10,-2.5) )
    target.sigma <- list( 1.5*diag(2),
                          2*matrix(c(2,0.5,0.5,1),2,2),
                          matrix(c(2.5,-1.5,-1.5,2.5),2,2) )
    Y <- rmvnorm(N/3,mean=target.mean[[1]],sigma=target.sigma[[1]])
    Y <- rbind(Y,rmvnorm(N/3,mean=target.mean[[2]],sigma=target.sigma[[2]]))
    Y <- rbind(Y,rmvnorm(N/3,mean=target.mean[[3]],sigma=target.sigma[[3]]))
    source.mean <- c(0,-5)
    source.sigma <- diag(x=1,nrow=d)
    X0 <- rmvnorm( n=N, mean=source.mean, sigma=source.sigma )
  } else {
    if( case == "Mx2G" ) {
      stopifnot( obj.fun %in% c("FisW2","FisSymKL") )
      Ff <- ifelse(obj.fun=="FisW2",F.sqdW2,F.symKL)
      set.seed(42)
      target.mean <- c(5,-5)
      target.sigma <- diag(2)
      Y <- rmvnorm( n=N, mean=target.mean, sigma=target.sigma )
      source.mean <- list( c(0,5), c(5,5), c(12,0) )
      source.sigma <- list( matrix(c(2,0.5,0.5,2),2,2),
                            2*matrix(c(2,-0.5,-0.5,1),2,2),
                            matrix(c(2,1,1,2),2,2) )
      X0 <- rmvnorm(N/3,mean=source.mean[[1]],sigma=source.sigma[[1]])
      X0 <- rbind(X0,rmvnorm(N/3,mean=source.mean[[2]],sigma=source.sigma[[2]]))
      X0 <- rbind(X0,rmvnorm(N/3,mean=source.mean[[3]],sigma=source.sigma[[3]]))
    } else {
      stop("ERROR!")
    }
  }
}



cat("> Starting constrained-JKO with eps =",eps,"\n")
elapsed <- Sys.time()
res <- cJKO( X=X0, Ff=Ff, eps=eps, K=K.cjko, patience=patience.cjko )
elapsed <- difftime(Sys.time(),elapsed,units="secs")
cat("\n> Solved in:",elapsed,"[secs]\n\n")

res$patience <- patience.cjko
res$eps <- eps
res$K <- K.cjko
res$X0 <- X0
res$Y <- Y
res$srcMeans <- source.mean
res$srcSigmas <- source.sigma
res$tgtMeans <- target.mean
res$tgtSigmas <- target.sigma

resFolder <- paste0("cJKO_results_N",N)
if( !dir.exists(resFolder) )
  dir.create(resFolder)
ix <- 1
while( file.exists( paste0(resFolder,"/cJKO_",case,"_",obj.fun,"_",ix,".RData") ) )
  ix <- ix+1
saveRDS( res, paste0(resFolder,"/cJKO_",case,"_",obj.fun,"_",ix,".RData") )
