## extra simulations ####
## multivariate normal
ud0 = c(0.5, 0.8, 1, 1.5)   ## mean vector for case
p = length(ud0)
un0 = rep(0, p)         ## mean vector for control
mu0 = ud0 - un0
cvmd0 = matrix(0.8, p, p)   ## covarince matrix for case
diag(cvmd0) = 1
cvmn0 = matrix(0.2, p, p)   ## covarince matrix for control
diag(cvmn0) = 1

result = max_pauc_simu(u = 0.1, n = 40, m = 40, type = 'mvn', 
                  method = 'ALG2', totalrun = 100)
mean(result$pauc)
sd(result$pauc)

roc(result$pool.case, result$pool.control, plot = TRUE)

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


## for real example
result = max_pauc_real(casedata, controldata, u = 0.01, method = 'ALG1', totalrun = 10)
mean(result$pauc)
roc1 = roc(cases = result$pool.case, controls = result$pool.control)

result = max_pauc_real(casedata, controldata, u = 0.01, method = 'LOGIT', totalrun = 10)
roc2 = roc(cases = result$pool.case, controls = result$pool.control)
mean(result$pauc)

## plot some ROC curves ####

plot(sort(1-roc1$specificities), sort(roc1$sensitivities), 
     xlim=c(0, 0.01), ylim=c(0, 0.8), type='l', xlab='1-Specificity', 
     ylab='Sensitivity', lwd=2, lty=1)
lines(sort(1-roc2$specificities), sort(roc2$sensitivities), 
      col=2, lwd=2, lty=2)


## use ggplot2
spe = sapply(rocs, function(x) sort(1-x$specificities))
FPR = matrix(spe, ncol=1, byrow=F)
sen = sapply(rocs, function(x) sort(x$sensitivities))
TPR = matrix(sen, ncol=1, byrow=F)
Methods = rep(mlist, each=nrow(spe))
datg = data.frame(FPR, TPR, Methods)
library(grid) # for unit
ggplot(datg, aes(x=FPR, y= TPR, group=Methods, 
                 colour=Methods, linetype = Methods)) + 
  xlab('1-Specificity') + ylab('Sensitivity') + 
  theme_bw() + theme(legend.position=c(.9,.25)) + geom_line(size=1.5) +
  xlim(c(0, 0.1)) + ylim(c(-0.05, 0.9)) + 
  theme(panel.grid = element_blank()) + theme(legend.key.size = unit(1.5, "cm")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 




## another version of implementation ####

## a simulation example
rept = 100   ## repeat 100 times
n = m = 40
u = 0.1      ## upper bound of fpr
type = 'mvn'            ## indicates multivariate normal

## setting the true mean vectors and covariance matrices
ud0 = c(0.5, 0.8, 1, 1.5)   ## mean vector for case
p = length(ud0)
un0 = rep(0, p)         ## mean vector for control
mu0 = ud0 - un0
cvmd0 = matrix(0.8, p, p)   ## covarince matrix for case
diag(cvmd0) = 1
cvmn0 = matrix(0.2, p, p)   ## covarince matrix for control
diag(cvmn0) = 1

pauc.train = pauc.test = rep(0, rept)
method = "LOGIT"       ## which method to obtain the linear combination
set.seed(1234)        ## for reproducability
for(i in 1:rept){
  ## generate data
  cases.train = data_simu(type, n, ud0, cvmd0)
  controls.train = data_simu(type, m, un0, cvmn0)
  
  cases.test = data_simu(type, n, ud0, cvmd0)
  controls.test = data_simu(type, m, un0, cvmn0)
  result = MaximzPAUC(cases.train, controls.train, u, method,
                         cases.test, controls.test)
  pauc.train[i] = result$pauc.train
  pauc.test[i] = result$pauc.test
}
#round(mean(pauc.train), 4)
#round(sd(pauc.train), 4)
round(mean(pauc.test), 4)
round(sd(pauc.test), 4)



