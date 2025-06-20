# Title        : Script for the MLE of lineage frequencies and MOI parameter using Incomplete Data Model
# Objective    : Estimate lineage frequencies, MOI,
# Created by   : Meraj Hashemi, Kristan A. Schneider
# Created on   : 25.04.22
# Last modified: 16.05.23


DatImp <- function(path){
  if(substring(path,nchar(path)-3,nchar(path))==".xls"){
    dat <- openxlsx::read.xlsx(path,1)
  }
  else{
    if(substring(path,nchar(path)-4,nchar(path))==".xlsx"){
      dat <- openxlsx::read.xlsx(path,1)
    }
    else{
      if(substring(path,nchar(path)-3,nchar(path))==".txt"){
        dat <- read.table(path,header=TRUE, sep="\t")
      }
      else{
        if(substring(path,nchar(path)-3,nchar(path))==".csv"){
          dat <- read.csv(path,header=TRUE,sep=";")
        }
      }
    }  
  }
  dat
}  

#************************************************************************************
#This function calculate Nk
#************************************************************************************    

Nk <- function(dat){
  for(k in 1:nrow(dat)){
    if(dat[k,1]==""||is.na(dat[k,1])){
      dat[k,1] <- dat[k-1,1]
    }
  }
  N <- length(unique(dat[,1]))
  dat <- dat[!is.na(dat[,2]),]
  Nplus <- length(unique(dat[,1]))
  Nknum <- length(unique(dat[,2]))
  dat <- dat[!duplicated(dat),]
  out <- list(N,t(as.matrix(table(dat[,2]))),N-Nplus)
  names(out) <- c("N","N_k","n_0")
  out
}

########################################################################################
#------------------------- Functions for calculating the MLEs---------------------------
########################################################################################

#' Funktion to derive MLE for the IDM or OM
#'
#' @param N integer-valued float; Sample size
#' @param n_0 integer-valued float; number of empty records
#' @param N_k float vector integer valued;each component corresponds to the number of
#'   times a lineage is found in the dataset
#' @param lambda_initial float; initial value of lambda for numerical algorithm, it should only be adjusted if there is prior 
#'  information on the the value of lambda or if numerical problems occur
#' @param eps_initial float; initial value of epsilon (probability of the
#'   lineages remain undetected) used by the numerical algorithm.  It should only be adjusted if there is prior 
#'  information on the the value of lambda or if numerical problems occur.
#'
#' @return the function returns a list of values as follows:
#'         1) the MLE of the probability of lineages remaining undetected (epsilon) - this output is omitted if option model="OM" is specified;
#'         2) the MLE of the MOI parameter (lambda);
#'         3) the MLE of the average MOI (psi);
#'         4) the MLE of the lineage frequencies;
#'         5) the inverse Fisher information (estimates for the parameter epsilon are omitted if option model="OM" is specified).
#'
#' @examples MLE_IDM(40, 1, c(23,27), 1, 0.1)
MLE <- function(N, N_k, n_0=0, model = "IDM", lambda_initial = 1, eps_initial=0.1){
  eps <- 1e-12
  n <- length(N_k)
  if(!is.numeric(N)){
    warning("Argument N must be an interge valued float")
  }else if(!is.numeric(N_k)){
    warning("Argument N_k must be an interge valued float vector")
  }else if(!is.numeric(n_0)){
    warning("Argument n_0 must be an interge valued float")
  }else{
    N1 <- floor(N)
    N_k1 <- floor(N_k)
    n_01 <- floor(n_0)
    if(N-N1>eps){
      warning(paste("Argument N must be an interge valued float, it was changed to N=", N1,sep=""))
      N <- N1
    }else if(sum(N_k-N_k1)>eps){
      warning(paste("2- Argument N_k must be an interge valued float, it was changed to N_k=", N_k1,sep=""))
      N_k <- N_k1
    }else if(n_0 - n_01 >eps){
      warning(paste("Argument n_0 must be an interge valued float, it was changed to n_0=", n_01,sep=""))
      n_0 <- n_01
    }
    if(n_0 <0 || min(N_k)< 0 || N<0 || min(N-N_k)<0 || n_0>min(N-N_k) || (sum(N_k) <N && is.element(model,c("OM"))) || (sum(N_k) + n_0 <N  )){
      warning("The data does not satisfiy the requirements, N, N_k, n_0 must be natural numbers, max(N_k)<N, and n_0 <= min(N-N_k), n_0+N1,+ .. + N_n >= N")
    }else{
      if(n==1){
        final <- list(n_0/N,NA,NA,1,NA,NA)
        names(final) <- c("probability of lineages remain undetected", "MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
        if(is.element(model,c("OM"))){
          final <- list(NA,NA,1,NA,NA)
          names(final) <- c("MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
        }
        final
      }else{
        if(is.element(model,c("OM"))){
          N <- N-n_0
          n_0 <-0
        }
        if(n_0==0){
          if(sum(N_k)==N){
            final <- list(0,0,1,N_k/N,NA,NA)
            names(final) <- c("probability of lineages remain undetected", "MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
            if(is.element(model,c("OM"))){
              final <- list(0,1,N_k/N,NA,NA)
              names(final) <- c( "MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
            }
            final
          }else if(min(N-N_k)==0){
            final <- list(0,Inf,Inf,NA,NA,NA)
            names(final) <- c("probability of lineages remain undetected", "MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
            if(is.element(model,c("OM"))){
              final <- list(Inf,Inf,NA,NA,NA)
              names(final) <- c("MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
            }
            final
          }else{
           
            final <- MLE1(N, N_k, n_0, model, lambda_initial, eps_initial)
            final
          }
        }else{ #n_0>0
          if(prod(1-N_k/N)<=n_0/N){
            final <- list(min(1-N_k/N),Inf,Inf,NA,NA,NA)
            names(final) <- c("probability of lineages remain undetected", "MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
            if(is.element(model,c("OM"))){
               final <- list(Inf,Inf,NA,NA,NA)
               names(final) <- c("MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
            }
            final
          }else if(prod(1-N_k/N)>n_0/N){
            if(sum(N_k)==N-n_0){
              final <- list(n_0/N,0,1,NA,NA,NA)
              names(final) <- c("probability of lineages remain undetected", "MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
#              if(is.element(model,c("OM"))){
#                final <- list(0,1,N_k/N,NA,NA)
#                names(final) <- c( "MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
#              }
              final
            }else{
               final <- MLE1(N, N_k, n_0, model, lambda_initial, eps_initial)
               final
            }
          }
        }
      }
    }
  }
}  


#---------------------------------internal function-------------------------------

MLE1 <- function(N, N_k, n_0=0, model = "IDM", lambda_initial = 1, eps_initial=0.1){
  if(model == "OM"){
    if(n_0>0){
      print(paste("Option model= `OM' neglects n_0=",n_0," samples and adjusts sample size to N=",N-n_0,sep=""))
    }
    MLE_OM(N-n_0,N_k,lambda_initial)
  }else if(model == "IDM"){
    if(n_0>0){
      MLE_IDM(N, N_k,n_0,lambda_initial, eps_initial)
    }else{
      inp <- MLE_OM(N-n_0,N_k,lambda_initial)
      FI <- inp[[4]]
      n <- length(inp[[3]])
      FInf <- array(NA,c(n+2,n+2))
      pick <- c(TRUE,FALSE,rep(TRUE,n))
      #print(pick)
      FInf[pick,pick] <- inp[[4]]
      nam <-c("lam","eps",paste("p",1:n,sep="."))
      colnames(FInf) <- nam
      rownames(FInf) <- nam
      
      lam <- inp[[2]]
      el <- exp(lam)
      adj <- el*(el-lam-1)/(el-1)^2
      FInfadj <- FInf
      FInfadj[1,] <- FInf[1,]*adj
      FInfadj[,1] <-  FInf[,1]*adj
      nam <-c("psi","eps",paste("p",1:n,sep="."))
      colnames(FInfadj) <- nam
      rownames(FInfadj) <- nam

      final <- list(0,inp[[1]],inp[[2]],inp[[3]],FInf,FInfadj)
      names(final) <- c("probability of lineages remain undetected", "MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
      final
    }
        
  }else{
    warning("option model needs to be eiter `IDM' or `OM'")
  }
}    

################################ The MLE of the IDM #################################

#' The EM algorithm to derive the MLE for the IDM
#'
#' @param N integer-valued float; Sample size
#' @param n0 integer-valued float; number of empty records
#' @param Nk float vector integer valued;each component corresponds to the number of
#'   times a lineage is found in the dataset
#' @param lambda_initial float; initial value of lambda
#' @param eps_initial float; initial value of epsilon (probability of the
#'   lineages remain undetected)
#'
#' @return the function returns a list of values as follows:
#'         1) the MLE of the probability of lineages remaining undetected (epsilon);
#'         2) the MLE of the MOI parameter (lambda);
#'         3) the MLE of the average MOI (psi);
#'         4) the MLE of the lineage frequencies;
#'         5) the inverse Fisher information evaluated at the LME.
#'
#' @examples MLE_IDM(40, 1, c(23,27), 1, 0.1)
MLE_IDM <- function(N, Nk,n0, lambda_initial, eps_initial) {
  thr1 <- 10^-8
  thr2 <- 10^-12
  thr3 <- 10^-20
  z <- Nk
  sel <- Nk!=0
  Nk <- Nk[sel]
  n <- length(Nk)
  Nnk <- N - Nk
  snk <- sum(Nk)
  #initial values
  lamt <- 0.1
  lamnext <- lambda_initial
  pkt <- as.vector(array(1/n,c(1,n))) + 0.1
  pnext <- as.vector(array(1/n,c(1,n)))
  epst <- 0.01
  epsnext <- eps_initial 
  while(abs(lamt - lamnext) +abs(epst - epsnext) + sqrt(sum((pkt - pnext)^2)) > thr1) {
    lamt <- lamnext
    epst <- epsnext
    pkt <- pnext
    nextiter <- EM_next_iteration(lamnext, pnext, epst, N, n0, Nk, Nnk, snk)
    pnext <- nextiter[[1]]
    epsnext<- nextiter[[2]] 
    wt <- nextiter[[3]]
    ntt <- wt/N
    lamnext <-  ntt + 1
    lamnext <- EM_lambda_Newton(lamnext, N, thr1, ntt)
  }
  
  pnextz <- array(0,length(z))  
  pnextz[z > 0] <- pnext 
  psi <- lamnext/(1-exp(-lamnext))
  pp <- array(0,length(sel))

  pick <-c( TRUE,epsnext>0, sel)
  FI <- FI(N,lamnext,pnext,epsnext)
  pick <-c( TRUE,TRUE, sel)
  n <- length(sel)
  FInf <- array(NA,c(n+2,n+2))
  FInf[pick,pick] <- FI
  nam <-c("lam","eps",paste("p",1:n,sep="."))
  colnames(FInf) <- nam
  rownames(FInf) <- nam

  
  
  lam <- lamnext
  el <- exp(lam)
  adj <- el*(el-lam-1)/(el-1)^2
  FInfadj <- FInf
  FInfadj[1,] <- FInf[1,]*adj
  FInfadj[,1] <-  FInf[,1]*adj
  nam <-c("psi","eps",paste("p",1:n,sep="."))
  colnames(FInfadj) <- nam
  rownames(FInfadj) <- nam

  final <- list(epsnext, lamnext, psi, pnextz,FInf,FInfadj)
  names(final) <- c("probability of lineages remain undetected", "MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
  final
}

#---------------------------------internal function-------------------------------

EM_lambda_Newton <- function(initial, N, thr, ntt) {
  lamt <- 0
  lamnext <- initial
  while (abs(lamt - lamnext) > thr) {
    lamt <- lamnext
    exp_l <- 1 - exp(-lamt)
    newt <- 1 - exp(-lamt)*ntt
    lamnext <- (ntt*(1 - exp(-lamt)*(1 + lamt)))/newt 
  }
  
  lamnext
}

#---------------------------------internal function-------------------------------


EM_next_iteration <- function (lamnext, pkt, epst, N, N0, Nk, Nnk, snk) {
  #prereuisite
  expt <- exp(lamnext*pkt)
  exp1t <- expt - 1
  exp1expt <- expt/exp1t
  expepst <- epst*exp1t + 1
  if (N0 > 0) {
    tt <- N0/(-1 + prod(expepst))
  }
  else {
    tt <- 0
  }
  
  exp2epst <- expt/expepst
  
  wt <- lamnext*( sum(pkt*(Nk*exp1expt + Nnk*epst*exp2epst)) 
                  + epst*tt*sum(pkt*exp2epst))
  exp1e <- exp1t/expepst
  vkt <- Nnk*exp1e
  vt <- epst*(sum(vkt) + sum(exp1e)*tt)
  ukt <- pkt*lamnext*(Nk*exp1expt + epst*exp2epst*(Nnk + tt))
  #next iteration
  pnext <- ukt/sum(ukt)
  epsnext <- 1/(1 + (snk/vt))
  list(pnext, epsnext, wt)
}


#---------------------------------internal function-------------------------------

#' The EM algorithm to derive the MLE for the IDM
#'
#' @param N integer-valued ; Sample size
#' @param lam float; MOI parameter
#' @param p float vector; vector of lineage frequencies
#' @param eps float; probability of the lineages remain undetected, default eps=0
#'
#' @return the inverse Fisher information matrix
#'
#' @examples FI(100, 1.1, c(0.5,0.45,0.05), 0.1)
#' 
FI <- function (N,lam,p,eps=0) {
  #prerequisite
  n <-length(p)
  if(eps==0){
    FI <- array(0,c(n+2,n+2))
    iFI <- array(0,c(n+2,n+2))
    el <- exp(lam)
    eml <-1/el
    elk <- exp(lam*p)
    eplkmo <- (elk-1)
    
    
    ## d^2L/dl^2
    FI[1,1] <- -1/(1-eml)+el*sum( p^2/eplkmo)
    
    ## d^2L/dldpk
    
    FI[1,1:n+2] <- - el*(1- p/eplkmo*lam )
    FI[1:n+2,1] <- FI[1,1:n+2]
    
    ## d^2L/dpk^2
    
    D <- (el*lam^2/eplkmo)
    #FI <- FI + diag(c(0,0,D ))
    
    FI <- N/(el-1)*FI
    FI[1:n+2,2] <- 1
    FI[2,1:n+2] <- FI[1:n+2,2] 
    
     ## use blockwise inverision
    D <- N/(el-1)*D
    D1 <- 1/D
    d <- -sum(D1)
    D2 <- c(-1,D1)
    D3 <- (D2%*%t(D2))/d + diag(c(0,D1))
    
    A0 <-FI[1,1]
    B <- FI[1,2:(n+2)]
    A1 <- FI[1:n+2,1:n+2]
    DC <- D3 %*% B
    A1 <- solve(A0 - (B %*% DC))
    
    iFI[1,1] <- A1
    C2 <-  -DC %*% A1 
    iFI[0:n+2,1] <- C2
    B2 <- t(C2)
    iFI[1,0:n+2] <- B2
    iFI[0:n+2,0:n+2] <- D3 - DC %*% B2 
    iFI <- iFI[-2,-2]
    
  }else{
    FI <- array(0,c(n+3,n+3))
    iFI <- array(0,c(n+3,n+3))
    el <- exp(lam)
    eml <-1/el
    elk <- exp(lam*p)
    Ak <- (elk-1)*eps +1
    tau <- prod(Ak)
    tau <- tau/(tau-1)
    pkAk <- t(p/Ak)
    eplkmo <- (elk-1)
    eplmoAl <- t(eplkmo/Ak)
    Tl <- eps * p * elk /Ak
    Te <- eplkmo/Ak
    Tp <- eps*lam*elk/Ak
    ## d^2L/dl^2
    FI[1,1] <- -1/(1-eml)+(1-eps)*el*(pkAk%*% (p/eplkmo)) + tau*sum(Tl)^2
    
    ## d^2L/dlde
    FI[1,2] <- -el*sum(pkAk)  + tau * sum(Tl)*sum(Te)
    FI[2,1] <- FI[1,2]
    
    ## d^2L/dldpk
    
    FI[1,1:n+3] <- - el*( 1- pkAk/eplkmo*lam*(1-eps) ) + tau * Tp *sum(Tl)
    FI[1:n+3,1] <- FI[1,1:n+3]
    
    ## d^2L/deps^2
    
    FI[2,2] <- el*sum(eplmoAl)/(1-eps)  + tau*sum(Te)^2
    
    ## d^2L/depsdpk
    
    FI[2,1:n+3] <- - el*lam/Ak  + tau*sum(Te) * Tp
    FI[1:n+3,2] <- FI[2,1:n+3]

    ## d^2L/dpk^2
    
    D <- (el*(1-eps)*lam^2/Ak/eplkmo + tau *  Tp^2)
    FI <- FI + diag(c(0,0,0,D ))
    
    ## use blockwise inverision
    FI <- N/(el-1)*FI
    
    D <- N/(el-1)*D
    D1 <- 1/D
    d <- -sum(D1)
    D2 <- c(-1,D1)
    D3 <- (D2%*%t(D2))/d + diag(c(0,D1))
    
    A0 <-FI[1:2,1:2]
    B <- FI[1:2,3:(n+3)]
    A1 <- FI[1:n+3,1:n+3]
    DC <- D3 %*% t(B)
    A1 <- solve(A0 - (B %*% DC))
    
    iFI[1:2,1:2] <- A1
    C2 <- - DC %*% A1 
    iFI[0:n+3,1:2] <- C2
    B2 <- t(C2)
    iFI[1:2,0:n+3] <- B2
    iFI[0:n+3,0:n+3] <- D3 - DC %*% B2 
    iFI <- iFI[-3,-3]
  }
  iFI
}

################################ The MLE of the OM #################################

#' function to derive the MLE for the original model
#'
#' @param N integer; Sample size
#' @param Nk vector of integers;each component corresponds to the number of
#'   times a lineage is found in the dataset
#' @param la float; initial value of lambda
#'
#' @return the function returns a list of values as follows: 1) the MLE of the
#'   MOI parameter (lambda) 2) the MLE of the average MOI (psi) 3) the MLE of
#'   the lineage frequencies
#' @export
#'
#' @examples MLE(97,c(22,25,49,32,18))
MLE_OM <- function(N,Nk,la=1){
  sel <- Nk
  Nk <- sel[sel>0]
  nk <- Nk/N
  l1 <- 2.5         # initial value
  l0 <- 0
  eps <- 10^(-8)       # precision 
  k <- 1
  while(abs(l0-l1)>eps && k<50 && l1>0){
    k <- k+1
    l0 <- l1
    l1 <- l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
  }
  if(k==50 || l1<0){
    for(st in 1:10){
      l1 <- st
      l0 <- l1+1
      k <- 1
      while(abs(l0-l1)>eps && k<100 && l1>0){
        k <- k+1
        l0 <- l1
        l1 <- l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
      }
      if(abs(l0-l1)<eps){
        break
      }
    }
    if(abs(l0-l1)>eps){
      l1 <- Rmpfr::mpfr(10*la,precBits=100)
      l0 <- l1+1
      while(abs(l0-l1)>eps){
        l0 <- l1
        l1=l0-(l0+sum(log(1-nk*(1-exp(-l0)))))/(1-sum(nk/(exp(l0)*(1-nk)+nk)))
      }
    }        
  }
  pk <- -1/l1*log(1-nk*(1-exp(-l1))) 
  n <- length(sel)
  pk1 <- array(0,n)
  pick <- sel>0
  pk1[pick] <- pk  
  psi <- l1/(1-exp(-l1))
  
  FI <- FI(N,l1,pk,0)
  pick <-c( TRUE, pick)
  FInf <- array(NA,c(n+1,n+1))
  FInf[pick,pick] <- FI
  nam <-c("lam",paste("p",1:n,sep="."))
  colnames(FInf) <- nam
  rownames(FInf) <- nam
  
  lam <- l1
  el <- exp(lam)
  adj <- el*(el-lam-1)/(el-1)^2
  FInfadj <- FInf
  FInfadj[1,] <- FInf[1,]*adj
  FInfadj[,1] <-  FInf[,1]*adj
  nam <-c("psi",paste("p",1:n,sep="."))
  colnames(FInfadj) <- nam
  rownames(FInfadj) <- nam

  
  out <- list(l1,psi,pk1,FInf,FInfadj)
  names(out) <- c("MOI parameter lambda","average MOI","lineage frequencies","inverse Fisher information","inverse Fisher information adjusted for average MOI")
  out	
}


########################################################################################
#------------------------- Functions for the simulation study --------------------------
########################################################################################

####################################### cpoiss #########################################

#' Generates conditional Poisson random numbers
#'
#' @param lambda float; the MOI parameter
#' @param N integer; the sample size 
#'
#' @return a vector of randomly generated conditional Poisson numbers
#'
#' @examples  cpoiss(1.5, 10)
#' 
cpoiss<-function(lambda,N){
  m <- 100 # to accelerate computation it is assumed that m<100 is generically drawn
  out <- rep(0,N)
  x <- runif(N,min=0,max=1)
  p0 <- ppois(0,lambda)
  nc <- 1/(1-exp(-lambda))
  pvec <- (ppois(1:m,lambda)-p0)*nc
  pvec <- c(pvec,1) 
  for (i in 1:N){
    k <- 1
    while(x[i] > pvec[k]){
      k <- k+1
    }
    if(k==m){ # if a m>=100 is drawn this is executed
      k <- k+1
      a <- dpois(k,lambda)*nc
      b <- pvec[m]+a
      while(x[i]>b){
        k <- k+1
        a <- a*lambda/k
        b <- b+a
      }
    }
    out[i] <- k
  }
  out
}

####################################### mnom #########################################

#' Generates molecular dataset for a given set of model parameters
#'
#' @param M either a positive integer or a vector of positive integers
#'   corresponding to conditional Poisson random numbers
#' @param p vector; vector of lineage frequencies
#'
#' @return a 0-1 matrix of size N x n where each row corresponds to a sample
#'   with
#'
#' @examples 
#'  mnom(8, c(0.25,0.25,0.25,0.25))
#'  
#'  mnom(c(8,5,6), c(0.25,0.25,0.25,0.25))
#' 
mnom <- function(M, p){
  N <- length(M)
  out<-matrix(0, N, length(p))
  for(k in 1:N){
    out[k,] <- rmultinom(1,M[k],p)
  }
  out <- out
  out
}


####################################### runsim #########################################

#' Generates a molecular dataset with incomplete information 
#'
#' @param data matrix; a 0-1 matrix corresponding to N blood samples 
#' @param eps float; the probability of lineages remaining undetected
#' @param N integer; sample size
#' @param n integer; number of lineages
#'
#' @return a list with the following values:
#'            1) number of empty records
#'            2) the vector of observed prevalences
#'            3) dataset with incomplete information
#'
#' @examples
IncompleteData <- function(data, eps){
  N <- nrow(data)
  n <- ncol(data)
  ran <- runif(N*n)
  ran <- matrix((ran > eps)*1, N, n)
  ran*data
}
