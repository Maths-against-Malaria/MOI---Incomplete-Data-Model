# Title        : Script for the MLE of lineage frequencies and MOI parameter using Incomplete Data Model
# Objective    : Estimate lineage frequencies, MOI,
# Created by   : Meraj Hashemi, Kristan A. Schneider
# Created on   : 25.04.22
# Last modified: 07.03.23


library('ggplot2')
library('stringr')

################################# This function imports data #################################

DatImp<-function(path){
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

################################ This function calculate Nk #################################

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

MLE <- function(N, N_k, n_0=0, model = "IDM", lambda_initial = 1, eps_initial=0.1){
  if(model == "OM"){
    print(n_0)
    MLE_OM(N-n_0,N_k,lambda_initial)
  }else if(model == "IDM"){
      MLE_IDM(N, N_k,n_0,lambda_initial, eps_initial)
  }else{
    Warning("option model needs to be eiter `IDM' or `OM'")
  }
}    

################################ The MLE of the IDM #################################

#' The EM algorithm to derive the MLE for the IDM
#'
#' @param N integer; Sample size
#' @param N0 integer; number of empty records
#' @param Nk vector of integers;each component corresponds to the number of
#'   times a lineage is found in the dataset
#' @param lambda_initial float; initial value of lambda
#' @param eps_initial float; inital value of epsilon (probability of the
#'   lineages remain undetected)
#'
#' @return the function returns a list of values as follows:
#'         1) the MLE of the probability of lineages remaining undetected (epsilon)
#'         2) the MLE of the MOI parameter (lambda)
#'         3) the MLE of the average MOI (psi)
#'         4) the MLE of the lineage frequencies
#'
#' @examples EM(40, 1, c(23,27), 1, 0.1)
MLE_IDM <- function(N, Nk,n0, lambda_initial, eps_initial) {
  thr1 <- 10^-8
  thr2 <- 10^-12
  thr3 <- 10^-20
  z <- Nk
  Nk <- Nk[Nk!=0]
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
  final <- list(epsnext, lamnext, psi, pnextz)
  #names(final) <- c("")
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
#' @examples MLE(c(22,25,49,32,18),97)
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
  pk1 <- array(0,length(sel))  
  pk1[sel>0] <- pk  
  psi <- l1/(1-exp(-l1))
  out <- list(l1,psi,pk1)
  names(out) <- c("MOI parameter lambda","average MOI","lineage frequencies")
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
#'
#' @return dataset with incomplete information
#'
#' @examples
IncompleteData <- function(data, eps){
  N <- nrow(data)
  n <- ncol(data)
  ran <- runif(N*n)
  ran <- matrix((ran > eps)*1, N, n)
  ran*data
}

