##*********************************************************##
####         part 2: two simulation examples             ####
##*********************************************************##

library(MaxmzpAUC)
library(copula)

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
method = "ALG1"       ## which method to obtain the linear combination
set.seed(1234)        ## for reproducability
for(i in 1:rept){
  ## generate data
  cases.train = data_simu(type, n, ud, cvmd)
  controls.train = data_simu(type, m, un, cvmn)
  
  cases.test = data_simu(type, n, ud, cvmd)
  controls.test = data_simu(type, m, un, cvmn)
  result = MaxmzpAUC(cases.train, controls.train, u, method,
                      cases.test, controls.test)
  pauc.train[i] = result$pauc.train
  pauc.test[i] = result$pauc.test
}

round(mean(pauc.test), 4)
round(sd(pauc.test), 4)

## MSE
nums = c(100, 150, 200, 400, 1000, 2000)
rept = 100
type = 'mvn'
len = length(nums)
err1 = err2 = matrix(0, rept, len)
for(i in 1:len){
  set.seed(1234)
  n = m = nums[i]
  for(j in 1:rept){
    cases.train = data_simu(type, n, ud, cvmd)
    controls.train = data_simu(type, m, un, cvmn)
    
    #cases.test = data_simu(type, n, ud, cvmd)
    #controls.test = data_simu(type, m, un, cvmn)
    result = MaxmzpAUC(cases.train, controls.train, u, 'ALG1')
    err1[j, i] = sum((tcoefs -  result$coefs)^2)
    
    result = MaxmzpAUC(cases.train, controls.train, u, 'ALG2')
    err2[j, i] =  sum((tcoefs -  result$coefs)^2)
  }
}

y1 = apply(err1, 2, function(x) sum(x^(1)))/rept
y2 = apply(err2, 2, function(x) sum(x^(1)))/rept
plot(nums[-1], y1[-1], lwd = 2, lty = 1, type = 'l', xlab = 'sample size', ylab = 'MSE')
lines(nums[-1], y2[-1], lwd = 2, lty = 2)
legend(x = 1600, y = 0.05, legend = c('ALG1', 'ALG2'), lty = c(1:2), 
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
controldata = log(controldata + 10)
casedata = log(casedata + 10)

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
rept = 100
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
    result = MaxmzpAUC(cases.train, controls.train, u, method = mlist[j],
                        cases.test, controls.test, greed = FALSE)
    pauc.train[i, j] = result$pauc.train
    pauc.test[i, j] = result$pauc.test
    #pool.case = c(pool.case, result$case.score)
    #pool.control = c(pool.control, result$control.score)
  }
  #rocs[[j]] = roc(cases = pool.case, controls = pool.control)
}

me_pauc = round(apply(pauc.test, 2, mean), 4)
sd_pauc = round(apply(pauc.test, 2, sd), 4)

dat = data.frame(mlist, me_pauc, sd_pauc)
dat


## tring a inverse normal transformation -- not work
dat = rbind(casedata, controldata)
dat = apply(dat, 2, function(x) qnorm((rank(x, ties = 'random')-0.1)/(n+m)))
casedata = dat[1:n, ]
controldata = dat[-(1:n), ]

## tring a box-cox transformation ####

outcome = as.factor(rep(c(1, 0), c(n, m)))

## box-cox transformation given lambda
fun <- function(x, lambda){
  ll = length(x)
  if(!all(x > 0)) x = x - min(x) + 0.5
  pd = (prod(x))^(1/ll)
  if(lambda==0){
    tx = pd*log(x)
  }else{
    tx = (x^lambda-1)/(lambda*pd^(lambda-1))
  }
  return(tx)
}


# find best lambda for box-cox trans for covariates
# in logistic model
box.cox.logit  <- function(lambdas, dat, outcome){
  
  
  len = length(lambdas)
  dev = rep(0, len)
  for(i in 1:len){
    tdat = apply(dat, 2, fun, lambdas[i])
    logitmod = glm(outcome ~., data = data.frame(tdat, outcome), 
                   family = binomial("logit"))
    dev[i] = logitmod$deviance
  }
  lbd = lambdas[which.min(dev)]
  tdat = apply(dat, 2, fun, lbd)
  
  return(list(tdat, lbd))
}

library(mvtnorm)
# for maximizing multivariate normal log-likelihood
box.cox.twog <- function(lambdas, cases, controls){
  len = length(lambdas)
  logl = rep(0, len)
  for(i in 1:len){
    tcases = apply(cases, 2, fun, lambdas[i] )
    tcontrols = apply(controls, 2, fun, lambdas[i])
    logl1 = logl2 = 1
    for(j in 1:nrow(tcases)){
      logl1 = logl1 * dmvnorm(tcases[j, ], mean = colMeans(tcases), 
                      sigma = cov(tcases))
    }
    for(j in 1:nrow(tcontrols)){
      logl2 = logl2 * dmvnorm(tcontrols[j, ], mean = colMeans(tcontrols), 
                      sigma = cov(tcontrols))
    }
    logl[i] = logl1 * logl2
  }
  lbd = lambdas[which.max(logl)]
  tcases = apply(cases, 2, fun, lbd )
  tcontrols = apply(controls, 2, fun, lbd)
  return(list(tcases, tcontrols, lbd))
}


# maximizing pooled likelihood
box.cox <- function(lambdas, dat){
  len = length(lambdas)
  logl = rep(1, len)
  for(i in 1:len){
    tdat = apply(dat, 2, fun, lambdas[i] )
    
    logl1 = 1
    for(j in 1:nrow(dat)){
      logl1 = logl1 * dmvnorm(tdat[j, ], mean = colMeans(tdat), 
                              sigma = cov(tdat))
    }
    
    logl[i] = logl1 
  }
  lbd = lambdas[which.max(logl)]
  tdat = apply(dat, 2, fun, lbd )
  return(list(tdat, lbd))
}


dat = rbind(casedata, controldata)
#dat = apply(dat, 2, function(x) x/max(abs(x)))
res = box.cox(seq(-3, 3, len = 300), dat)
res[[2]]
dat = res[[1]]
casedata = dat[1:n, ]
controldata = dat[-(1:n), ]
