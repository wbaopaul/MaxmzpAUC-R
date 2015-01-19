
##===============================================================##
## *****This program was wrote by Wenbao Yu for the paper********** 
## *****sumbmited to journal CSDA, titled:                   ******
## *****'Two simple algorithms on linear combination of multiple***
## ***** biomarkers to maximize partial area under the ROC curve'**
##===============================================================##


##*********************************************************##
####           part 1: list of all functions             ####
##*********************************************************##

library(MASS)
library(pROC)
library(copula)
library(ggplot2)

## function to get pauc by integration for a given combination
## for binormal ROC
pauc_integral <- function(u, alpha, mu, cvmn, cvmd){
  ua = sum(alpha * mu)/sqrt(t(alpha) %*% cvmd %*% alpha)
  va = sqrt(t(alpha) %*% cvmn %*% alpha)/sqrt(t(alpha) %*% cvmd %*% alpha)
  senf <- function(x){
    senc = pnorm(ua - va * qnorm(1 - x))
    return(senc)
  } 
  pauc = integrate(senf, 0, u)$value
  return(pauc)
}

## algorithm 1--semiparametric grid search
semigrid_pauc <- function(u, mu, cvmn, cvmd){
  ii = max(abs(cvmn - cvmd))
  if(ii <= 0.01) return(ginv(cvmn) %*% mu)
  
  b0 = seq(-1, 1, length = 100)
  len = length(b0)
  pauc = matrix(0, len, 4)
  for(i in 1:len){
    a1 = ginv(cvmd + b0[i]*cvmn) %*% mu
    a2 = ginv(-cvmd + b0[i]*cvmn) %*% mu
    a3 = ginv(b0[i]*cvmd + cvmn) %*% mu
    a4 = ginv(b0[i]*cvmd - cvmn) %*% mu
    pauc[i, 1] = pauc_integral(u, a1, mu, cvmn, cvmd)
    pauc[i, 2] = pauc_integral(u, a2, mu, cvmn, cvmd)
    pauc[i, 3] = pauc_integral(u, a3, mu, cvmn, cvmd)
    pauc[i, 4] = pauc_integral(u, a4, mu, cvmn, cvmd)
  }
  #pauc = round(pauc, 5)
  ind = which(pauc == max(pauc), arr.ind = T)[1, ]
  
  i = ind[1]
  j = ind[2]
  if(j == 1) w = c(1, b0[i])
  if(j == 2) w = c(-1, b0[i])
  if(j == 3) w = c(b0[i], 1)
  if(j == 4) w = c(b0[i], -1)
  
  coefs = ginv(w[1] * cvmd + w[2] *cvmn) %*% mu
  return(coefs)
  
}

## algorithm 2--parametric iterative method
iter_pauc <- function(a0, u, mu, cvmn, cvmd){
  run = 0
  b0 = a0
  repeat{
    run = run + 1
    qx = sum(t(a0) %*% cvmn %*% a0)
    qy = sum(t(a0) %*% cvmd %*% a0)
    v = sum(t(a0) %*% mu * sqrt(qx)/(qx + qy))
    sigma = sqrt(qy/(qx+qy))
    ct = qnorm(1 - u)
    c1 = sqrt(2 * pi) * sigma * pnorm((v - ct)/sigma)
    c2 = sigma^2 * exp(-(ct - v)^2/(2 * sigma^2))/(sqrt(qx) * qy)
    w1 = c1 * v/sqrt(qx) + c2 * qy
    w2 = c1 * v/sqrt(qx) - c2 * qx
    if (abs(w1) < 0.00001 && abs(w2) < 0.00001) return(b0)
    aa= ginv(w1 * cvmn + w2 * cvmd) %*% mu
    aa = aa/max(abs(aa))
    if(sum(abs(a0-aa)) < 0.01 || run == 100) break
    a0 = aa
  }
  return(a0)
}

## try p+1 different initial values for algorithm 2
multi_initial <- function(p, u, mu, cvmn, cvmd){
  a1 = ginv(cvmd + cvmn) %*% mu ## add Su-Liu's combination as an extra initial value
  ii = max(abs(cvmn - cvmd))
  if(ii <= 0.01) return(a1)
  
  a0 = rep(0, p)
  pauc = rep(0, p + 1)
  ## use e_i as the initial value
  for(i in 1:p){
    b0 = a0
    b0[i] = 1
    aa = iter_pauc(b0, u, mu, cvmn, cvmd) 
    pauc[i] = pauc_integral(u, aa, mu, cvmn, cvmd)
  }
  
  ## use su-liu's combination as initial value
  aan =  iter_pauc(a1, u, mu, cvmn, cvmd) 
  pauc[p+1] = pauc_integral(u, aan, mu, cvmn, cvmd)
  
  ind = which.max(pauc)
  if(ind != (p + 1)){
    a0[ind] = 1
    aa =  iter_pauc(a0, u, mu, cvmn, cvmd)
  }else{
    aa = aan 
  }
  return(aa)
}

## try many more different initial values -- optional
multi_initial_greed <- function(p, u, mu, cvmn, cvmd){
  a1 = ginv(cvmd + cvmn) %*% mu ## add Su-Liu's combination as an extra initial value
  ii = max(abs(cvmn - cvmd))
  if(ii <= 0.01) return(a1)
  ninits = 1000
  a0 = rep(0, ninits)
  pauc = rep(0, ninits + 1)
  ## use e_i as the initial value
  b0 = matrix(runif(p * ninits, -1, 1), p, ninits)
  for(i in 1:ninits){
    aa = iter_pauc(b0[, i], u, mu, cvmn, cvmd) 
    pauc[i] = pauc_integral(u, aa, mu, cvmn, cvmd)
  }
  
  ## use su-liu's combination as initial value
  aan =  iter_pauc(a1, u, mu, cvmn, cvmd) 
  pauc[ninits + 1] = pauc_integral(u, aan, mu, cvmn, cvmd)
  
  ind = which.max(pauc)
  if(ind != (ninits + 1)){
    aa =  iter_pauc(b0[, ind], u, mu, cvmn, cvmd) 
  }else{
    aa = aan 
  }
  return(aa)
}

## generate data from correlated multivariate non-normal distributions
data_mixture <- function(n, m, r = 0.5){
  myCop.norm <- ellipCopula(family = "normal", dim = 4, dispstr = "ex", param = r)
  myMvd <- mvdc(copula = myCop.norm, margins = c("norm", "norm", 'exp',"gamma"), 
                paramMargins = list(list(mean = par1, sd = 1), list(mean = par2, sd = 1.5), 
                                    list(rate = par3), list(shape = par4, scale = 1)))
  
  cases = (rMvdc(n, myMvd))
  
  myCop.norm <- ellipCopula(family = "normal", dim = 4, dispstr = "ex", param = 1-r)
  myMvd <- mvdc(copula = myCop.norm, margins = c("norm", "norm", 'exp',"gamma"), 
                paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1), 
                                    list(rate = 1), list(shape=1, scale=1)))
  
  controls = (rMvdc(m, myMvd))                               
  
  
  tlist = list('cases' = cases, 'controls' = controls)
  return(tlist)
}

## generate data for simulation
data_simu <- function(type = 'mvn', ...){
  if(type == 'mvn'){
    dat = mvrnorm(...)
  }
  if(type == 'log-mvn'){
    dat = exp(mvrnorm(...))
  }
  if(type == 'mixture'){
    dat = data_mixture(...)
  }
  return(dat)
}

## maximize pauc by a given method and data
MaximzPAUC <- function(cases.train, controls.train, u = 0.1, method = 'ALG1',
                       cases.test = NULL, controls.test = NULL, greed = FALSE){
  
  p = ncol(cases.train)
  n = nrow(cases.train)
  m = nrow(controls.train)
  
  # standandization
  s0 = apply(rbind(controls.train, controls.test), 2, sd)
  u0 = apply(rbind(controls.train, controls.test), 2, mean)
  controls.train = scale(controls.train, center = u0, scale = s0)
  cases.train = scale(cases.train, center = u0, scale = s0)
  
  # compute sample covanriance and mean ####
  cvmd = cov(cases.train)
  cvmn = cov(controls.train)
  
  ud = apply(cases.train, 2, mean)
  un = apply(controls.train, 2, mean)
  mu = ud - un
  
  # implement several linear combinations and ####
  # compare their empirical ROC curves and pAUCs ####
  
  if(method == 'S-L'){
    # Su-Liu(1993)'s combination
    coefs = ginv(cvmd + cvmn) %*% mu
  }
  if(method == 'LIU'){
    # liu(2005)'s combination
    bb = eigen(cvmd)
    btemp = bb$vectors %*% diag(1/sqrt(bb$values)) %*% t(bb$vectors)
    cvm = (btemp %*% cvmn %*% btemp)
    cc = eigen(cvm)
    min_eigenvec = cc$vectors[, which.min(cc$values)]
    coefs = btemp %*% min_eigenvec
  }
  if(method == 'LOGIT'){
    # logisti regression
    outcome = as.factor(rep(c(1, 0), c(n, m)))
    logitmod = glm(outcome ~., data = data.frame(rbind(cases.train, controls.train), outcome), 
                   family = binomial("logit"))
    coefs = as.numeric(coef(logitmod)[-1])
  }
  if(method == 'ALG1'){
    # ALG1 -- semiparametric grid search
    coefs = semigrid_pauc(u, mu, cvmn, cvmd)
  }
  if(method == 'ALG2'){
    # ALG2 -- parametric iterative approach
    if(greed == TRUE){
      coefs = multi_initial_greed(p, u, mu, cvmn, cvmd)
    }else{
      coefs = multi_initial(p, u, mu, cvmn, cvmd)
    }
  }
  
  coefs = coefs/max(abs(coefs))
  
  case_score = cases.train %*% coefs
  control_score = controls.train %*% coefs
  pauc.train = roc(cases = case_score, controls = control_score,
                   partial.auc = c(1, 1 - u))$auc
  
  
  # standandization and apply to test data
  pauc.test = NULL
  if(length(cases.test) >0 && length(controls.test) >0 ){
    controls.test = scale(controls.test, center = u0, scale = s0)
    cases.test = scale(cases.test, center = u0, scale = s0)
    case_score = cases.test %*% coefs
    control_score = controls.test %*% coefs
    pauc.test = roc(cases = case_score, controls = control_score,
                    partial.auc = c(1, 1 - u))$auc
  }
  
  return(list('coefs' = coefs, 'pauc.train' = pauc.train,
              'pauc.test' = pauc.test))
}


##*********************************************************##
####         part 2: two simulation examples             ####
##*********************************************************##

## 1. multivariate normal or log-maltibariate normal
rept = 100   ## repeat times
n = m = 40
u = 0.1      ## upper bound of fpr
type = 'mvn'            ## indicates multivariate normal
                        ## 'log-mvn' indicates log-multivariate normal

## setting the true mean vectors and covariance matrices
ud = c(0.5, 0.8, 1, 1.5)   ## mean vector for case
p = length(ud)
un = rep(0, p)         ## mean vector for control
mu = ud - un
cvmd = matrix(0.8, p, p)   ## covarince matrix for case
diag(cvmd) = 1
cvmn = matrix(0.2, p, p)   ## covarince matrix for control
diag(cvmn) = 1

pauc.train = pauc.test = rep(0, rept)
method = "LOGIT"       ## which method to obtain the linear combination
set.seed(1234)        ## for reproducability
for(i in 1:rept){
  ## generate data
  cases.train = data_simu(type, n, ud, cvmd)
  controls.train = data_simu(type, m, un, cvmn)
  
  cases.test = data_simu(type, n, ud, cvmd)
  controls.test = data_simu(type, m, un, cvmn)
  result = MaximzPAUC(cases.train, controls.train, u, method,
                      cases.test, controls.test)
  pauc.train[i] = result$pauc.train
  pauc.test[i] = result$pauc.test
}

round(mean(pauc.test), 4)
round(sd(pauc.test), 4)


## 2. simulation for mixture distribution ####
# set distribution parameters (see data_mixture for detail)
par1 = 1
par2 = 1
par3 = 1/3
par4 = 1.5

pauc.train = pauc.test = rep(0, rept)
method = "ALG1"  
type = 'mixture'
set.seed(1234)      
for(i in 1:rept){
  ## generate data
  temp = data_simu(type, n, m, r = 0.8)
  cases.train = temp$cases
  controls.train = temp$controls
  
  temp = data_simu(type, n, m, r = 0.8)
  cases.test = temp$cases
  controls.test = temp$controls
  result = MaximzPAUC(cases.train, controls.train, u, method,
                      cases.test, controls.test)
  pauc.train[i] = result$pauc.train
  pauc.test[i] = result$pauc.test
}
#round(mean(pauc.train), 4)
#round(sd(pauc.train), 4)
round(mean(pauc.test), 4)
round(sd(pauc.test), 4)
