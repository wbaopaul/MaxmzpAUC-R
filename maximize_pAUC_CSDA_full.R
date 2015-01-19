
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

## 3-dimension exhaust search for the true 
exhaust_true <- function(u, mu, cvmn, cvmd){
  b0 = seq(-1, 1, by = 0.001)
  len = length(b0)
  llen =len * len
  aa = matrix(0, 6 * llen, 3)
  pauc = rep(0, 6 * llen)
  for(j in 1:len){
    for(k in 1:len){
      aa[k + len*(j-1), ] = c(1, b0[j], b0[k])
      aa[k + len*(j-1)+llen, ] = c(-1, b0[j], b0[k])
      aa[k + len*(j-1)+2*llen, ] = c(b0[j], 1, b0[k])
      aa[k + len*(j-1)+3*llen, ] = c(b0[j], -1, b0[k])
      aa[k + len*(j-1)+4*llen, ] = c(b0[j], b0[k], 1)
      aa[k + len*(j-1)+5*llen, ] = c(b0[j], b0[k], -1)
    }
  }
  for(i in 1:(6 * llen)){
    #sencs = ComputeSenc(aa[i, ], mu0, cvmn0, cvmd0)
    #pauc[i] =  Compute_trapAUC(fpr, sencs, u)
    pauc[i] = pauc_integral(u, aa[i, ], mu, cvmn, cvmd)
  }
  
  ind = which(pauc == max(pauc))
  
  return(aa[ind, ])
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
  w = rep(0, 2)
  repeat{
    run = run + 1
    qx = sum(t(a0) %*% cvmn %*% a0)
    qy = sum(t(a0) %*% cvmd %*% a0)
    v = sum(t(a0) %*% mu*sqrt(qx)/(qx+qy))
    sigma = sqrt(qy/(qx+qy))
    ct = qnorm(1 - u)
    c1 = sqrt(2 * pi) * sigma * pnorm((v - ct)/sigma)
    c2 = sigma^2 * exp(-(ct-v)^2/(2 * sigma^2))/(sqrt(qx) * qy)
    w[1] = c1 * v/sqrt(qx) + c2 * qy
    w[2] = c1 * v/sqrt(qx) - c2 * qx
    #if (max(abs(w)) < 0.001) return( ginv(cvmn+cvmd) %*% mu)
    if (max(abs(w)) < 0.001) return( b0)
    aa= ginv(w[1] * cvmn + w[2] * cvmd) %*% mu
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
  
  #a1 = semigrid_pauc(u, mu, cvmn, cvmd)
  a0 = rep(0, p)
  #a0 = a1
  pauc = rep(0, p + 1)
  ## use e_i as the initial value
  for(i in 1:p){
    b0 = a0
    b0[i] = 1
    #b0[i] = ifelse(mu[i] > 0, 1, -1)
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

## try much more random initial values -- optional
multi_initial_greed1 <- function(p, u, mu, cvmn, cvmd){
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

## tr more specific initial values
multi_initial_greed <- function(p, u, mu, cvmn, cvmd){
  asl = ginv(cvmd + cvmn) %*% mu ## add Su-Liu's combination as an extra initial value
  ii = max(abs(cvmn - cvmd))
  if(ii <= 0.01) return(asl)
  
  a1 = semigrid_pauc(u, mu, cvmn, cvmd)
  
  ninits = 2*p
  pauc = rep(0, ninits + 1)
  ## use e_i as the initial value
  b0 = matrix(0, p, ninits)
  for(i in 1:p){
    b0[i, i] = 1
  }
  for(i in (p+1):(2*p)){
    b0[, i] = a1
    b0[i-p, i] = 0
  }
  for(i in 1:ninits){
    aa = iter_pauc(b0[, i], u, mu, cvmn, cvmd) 
    pauc[i] = pauc_integral(u, aa, mu, cvmn, cvmd)
  }
  
  ## use su-liu's combination as initial value
  aan =  iter_pauc(asl, u, mu, cvmn, cvmd) 
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
MaxmzpAUC <- function(cases.train, controls.train, u = 0.1, method = 'ALG1',
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
  
  case.score = cases.train %*% coefs
  control.score = controls.train %*% coefs
  pauc.train = roc(cases = case.score, controls = control.score,
                   partial.auc = c(1, 1 - u))$auc
  
  
  # standandization and apply to test data
  pauc.test = NULL
  if(length(cases.test) >0 && length(controls.test) >0 ){
    controls.test = scale(controls.test, center = u0, scale = s0)
    cases.test = scale(cases.test, center = u0, scale = s0)
    case.score = cases.test %*% coefs
    control.score = controls.test %*% coefs
    pauc.test = roc(cases = case.score, controls = control.score,
                    partial.auc = c(1, 1 - u))$auc
  }
  
  return(list('coefs' = coefs, 'pauc.train' = pauc.train,
              'pauc.test' = pauc.test, 'case.score' = case.score, 
              'control.score' = control.score))
}


## some earlier functions ####
## for simulation -- maxmize pAUC for a given method
## the output is a list includes:
## the predictive pauc value and coefs for each run, the pooled scores for case and control
max_pauc_simu <- function(u, n, m, type = 'mvn', method = 'ALG1',  totalrun = 100, r = 0.8){
  
  set.seed(1234)
  pauc = rep(0, totalrun)
  pool.case = pool.control = NULL
  
  #p = length(mu0)
  p = 4
  allcoefs = matrix(0, p, totalrun)
  
  ## repeat totalrun times
  run = 0
  repeat{
    run = run + 1
    
    # get training data for each run
    if(type == 'mixture'){
      temp = data_simu(type, n, m, r)
      cases.train = temp$cases
      controls.train = temp$controls
      
      temp = data_simu(type, n, m, r)
      cases.test = temp$cases
      controls.test = temp$controls
    }else{
      cases.train = data_simu(type, n, ud0, cvmd0)
      controls.train = data_simu(type, m, un0, cvmn0)
      
      cases.test = data_simu(type, n, ud0, cvmd0)
      controls.test = data_simu(type, m, un0, cvmn0)
    }
    
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
    
    # standandization
    controls.test = scale(controls.test, center = u0, scale = s0)
    cases.test = scale(cases.test, center = u0, scale = s0)
    
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
      coefs = multi_initial(length(mu), u, mu, cvmn, cvmd)
    }
    
    coefs = coefs/max(abs(coefs))
    
    case_score = cases.test %*% coefs
    control_score = controls.test %*% coefs
    pauc[run] = roc(cases = case_score, controls = control_score,
                    partial.auc = c(1, 1 - u))$auc
    pool.case = c(pool.case, case_score)
    pool.control = c(pool.control, control_score)
    
    allcoefs[, run] = coefs
    
    if (run == totalrun) break
  }
  
  return(list('pauc' = pauc, 'allcoefs' = allcoefs,
              'pool.case' = pool.case, 'pool.control' = pool.control))
}

## for real examples
## cases, controls is the whole dataset for real example
max_pauc_real <- function(cases, controls, u, method = 'ALG1',  totalrun = 100, r = 0.8){
  p = ncol(cases)
  n = nrow(cases)
  m = nrow(controls)
  n1 = round(n/2)
  m1 = round(m/2)
  #fpr <<- seq(0, u, length = 100) ## generate enough fpr points for calculating pauc
  
  pauc = rep(0, totalrun)
  set.seed(1234)
  pool.case = pool.control = NULL
  allcoefs = matrix(0, p, totalrun)
  
  # standandization
  s0 = apply(controls, 2, sd)
  u0 = apply(controls, 2, mean)
  controls = scale(controls, center = u0, scale = s0)
  cases = scale(cases, center = u0, scale = s0)
  
  ## repeat totalrun times
  run = 0
  repeat{
    run = run + 1
    
    # get training and testing data for each run
    set1 = sample(1:n, n1)
    set2 = sample(1:m, m1)
    cases.train = cases[set1, ]
    controls.train = controls[set2, ]
    cases.test = cases[-set1, ]
    contol.test = controls[-set2, ]
    
    cases.test = as.matrix(cases.test)
    controls.test = as.matrix(contol.test)
    
    # compute sample covanriance and mean ####
    cvmd = cov(cases.train)
    cvmn = cov(controls.train)
    ud = apply(cases.train, 2, mean)
    un = apply(controls.train, 2, mean)
    mu = ud-un
    
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
      outcome = as.factor(rep(c(1, 0), c(n1, m1)))
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
      coefs = multi_initial(p, u, mu, cvmn, cvmd)
    }
    
    coefs = coefs/max(abs(coefs))
    
    case_score = cases.test %*% coefs
    control_score = controls.test %*% coefs
    pauc[run] = roc(cases = case_score, controls = control_score,
                    partial.auc = c(1, 1 - u))$auc
    pool.case = c(pool.case, case_score)
    pool.control = c(pool.control, control_score)
    
    allcoefs[, run] = coefs
    
    if (run == totalrun) break
  }
  
  return(list('pauc' = pauc, 'allcoefs' = allcoefs,
              'pool.case' = pool.case, 'pool.control' = pool.control))
}



##*********************************************************##
####         part 2: some examples             ####
##*********************************************************##


## 1. multivariate normal or log-maltibariate normal ####
rept = 100   ## repeat times
n = m = 40
u = 0.1      ## upper bound of fpr
type = 'mvn'            ## indicates multivariate normal
## 'log-mvn' indicates log-multivariate normal

## setting the true mean vectors and covariance matrices
ud = c(0.5, 0.8, 1)   ## mean vector for case
p = length(ud)
un = rep(0, p)         ## mean vector for control
mu = ud - un
cvmd = matrix(0.8, p, p)   ## covarince matrix for case
diag(cvmd) = 1
cvmn = matrix(0.2, p, p)   ## covarince matrix for control
diag(cvmn) = 1

pauc.train = pauc.test = rep(0, rept)
method = "ALG1"       ## which method to obtain the linear combination
set.seed(1234)        ## for reproducability
for(i in 1:rept){
  ## generate data
  cases.train = data_simu(type, n, ud, cvmd)
  controls.train = data_simu(type, m, un, cvmn)
  
  cases.test = data_simu(type, n, ud, cvmd)
  controls.test = data_simu(type, m, un, cvmn)
  result = MaximzpAUC(cases.train, controls.train, u, method,
                      cases.test, controls.test)
  pauc.train[i] = result$pauc.train
  pauc.test[i] = result$pauc.test
}
#round(mean(pauc.train), 4)
#round(sd(pauc.train), 4)
round(mean(pauc.test), 4)
round(sd(pauc.test), 4)


## run all method at once ####
mlist = c('S-L', 'LIU', 'LOGIT', 'ALG1', 'ALG2')
u = 0.1
n = m = 40
type = 'mvn'
me_pauc = sd_pauc = NULL
rocs  = list()
for(i in 1:5){
  result = max_pauc_simu(u, n, m , type,
                         mlist[i], totalrun = 100)
  me_pauc = c(me_pauc, mean(result$pauc))
  sd_pauc = c(sd_pauc, sd(result$pauc))
  rocs[[i]] = roc(cases = result$pool.case, controls = result$pool.control)
  
}
dat = data.frame(me_pauc, sd_pauc)
dat = round(dat, 4)
dat

write.table(paste('n=m=', n, 'u =', u, 'type=', type, 'mu1=', ud0[1], 'rho=', rho1),
            file = 'summary_rerun.txt', append = T)
write.table(dat, file = 'summary_rerun.txt', append = T, row.names = mlist)


## plot ROC curves
spe = sapply(rocs, function(x) sort(1 - x$specificities))
FPR = matrix(spe, ncol = 1, byrow = F)
sen = sapply(rocs, function(x) sort(x$sensitivities))
TPR = matrix(sen, ncol = 1, byrow = F)
Methods = rep(mlist, each=nrow(spe))
datg = data.frame(FPR, TPR, Methods)
library(grid) # for unit
ggplot(datg, aes(x = FPR, y = TPR, group = Methods, 
                 colour = Methods, linetype = Methods)) + 
  xlab('1-Specificity') + ylab('Sensitivity') + 
  theme_bw() + theme(legend.position =c(.9, .31)) + geom_line(size = 1.5) +
  xlim(c(0, 0.1)) + ylim(c(0, 0.85)) + 
  theme(panel.grid = element_blank()) + theme(legend.key.size = unit(1.5, "cm")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 



## search the true combination --- only for 3 markers
tcoefs = exhaust_true(u = 0.1, mu0, cvmn0, cvmd0)
tcoefs

## get l2 error for ALG1 and ALG2 ####
nums = c(100, 200, 500, 1000, 2000)
err1 = err2 = NULL
for(i in nums){
  n = m = i
  result = max_pauc_simu(u = 0.1, n, m, type = 'mvn', 
                         method = 'ALG1', totalrun = 100) 
  err1 = c(err1, median(apply(result$allcoefs, 2, function(x) sqrt(sum((tcoefs - x)^2)))))
  
  result = max_pauc_simu(u = 0.1, n, m, type = 'mvn', 
                         method = 'ALG2', totalrun = 100) 
  err2 = c(err2, median(apply(result$allcoefs, 2, function(x) sqrt(sum((tcoefs - x)^2)))))
  
}
err1
err2
plot(nums, err1, lwd = 2, lty = 1, type = 'l', xlab = 'sample size', ylab = 'L2 error')
lines(nums, err2, lwd = 2, lty = 2)
legend(x = 1600, y = 0.2, legend = c('ALG1', 'ALG2'), lty = c(1:2), 
       lwd = c(2, 2), bty = 'n')


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
  result = MaximzpAUC(cases.train, controls.train, u, method,
                      cases.test, controls.test)
  pauc.train[i] = result$pauc.train
  pauc.test[i] = result$pauc.test
}
#round(mean(pauc.train), 4)
#round(sd(pauc.train), 4)
round(mean(pauc.test), 4)
round(sd(pauc.test), 4)


## real examples ####
# set1
dat = read.csv('BreastTissue.csv')
nlabel = c('con', 'adi', 'gla')
controldata = dat[dat$Class%in%nlabel, -c(1:2)]
casedata = dat[!dat$Class%in%nlabel, -c(1:2)]
controldata = log(controldata+10)
casedata = log(casedata+10)

# set2 - magic gamma telescope data  -- good example
dat = read.table('magic04.txt', sep=',')
casedata = dat[dat[, 11] == 'g', -11]
controldata = dat[dat[, 11] != 'g', -11]

# set 3 - DMD ---- OK
casedata = read.table('DMD_carrier.txt')
controldata = read.table('DMD_normal.txt')
casedata = log(casedata[, 5:8])
controldata = log(controldata[, 5:8])
casedata = casedata[complete.cases(casedata), ]
controldata = controldata[complete.cases(controldata), ]

casedata = as.matrix(casedata)
controldata = as.matrix(controldata)

## run a real example ####
u = 0.01
rept = 10
n = nrow(casedata)
m = nrow(controldata)
n1 = floor(n/2)
m1 = floor(m/2)
mlist = c('S-L', 'LIU', 'LOGIT', 'ALG1', 'ALG2')
pauc.train = pauc.test = matrix(0, rept, length(mlist))
rocs = list()
for(j in 1:length(mlist)){
  set.seed(1234)        ## for reproducability
  pool.case = pool.control = NULL
  for(i in 1:rept){
    ## generate data
    set1 = sample(1:n, n1)
    set2 = sample(1:m, m1)
    cases.train = casedata[set1, ]
    controls.train = controldata[set2, ]
    
    cases.test = casedata[-set1, ]
    controls.test = controldata[-set2, ]
    result = MaximzpAUC(cases.train, controls.train, u, method = mlist[j],
                        cases.test, controls.test, greed = FALSE)
    pauc.train[i, j] = result$pauc.train
    pauc.test[i, j] = result$pauc.test
    pool.case = c(pool.case, result$case.score)
    pool.control = c(pool.control, result$control.score)
  }
  rocs[[j]] = roc(cases = pool.case, controls = pool.control)
}

me_pauc = round(apply(pauc.test, 2, mean), 4)
sd_pauc = round(apply(pauc.test, 2, sd), 4)

dat = data.frame(mlist, me_pauc, sd_pauc)
dat

## plot ROC curves
spe = sapply(rocs, function(x) sort(1 - x$specificities))
FPR = matrix(spe, ncol = 1, byrow = F)
sen = sapply(rocs, function(x) sort(x$sensitivities))
TPR = matrix(sen, ncol = 1, byrow = F) 
Methods = rep(mlist, each=nrow(spe))
datg = data.frame(FPR, TPR, Methods)
library(grid) # for unit
ggplot(datg, aes(x = FPR, y = TPR, group = Methods, 
                 colour = Methods, linetype = Methods)) + 
  xlab('1-Specificity') + ylab('Sensitivity') + 
  theme_bw() + theme(legend.position =c(.9, .27)) + geom_line(size = 1.5) +
  xlim(c(0, 0.1)) + ylim(c(-0.14, 0.7)) + 
  theme(panel.grid = element_blank()) + theme(legend.key.size = unit(1.5, "cm")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
