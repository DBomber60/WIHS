rm(list=ls())
source('drivers.logistic.R')

dat2 = read.csv("combined2.csv")[,-1]
dat2$VL_0 = ifelse(dat2$VLOAD_0 < exp(10), 0, 1)
dat2$VL_1 = ifelse(dat2$VLOAD_1 < exp(10), 0, 1)

# logpi(beta =c(1,1,1), y = dat2$VL_1, X = cbind(1, dat2$A2_0, dat2$VL_0), dstar = c(1,1,1))
# grad(beta =c(1,1,1), y = dat2$VL_1, X = cbind(1, dat2$A2_0, dat2$VL_0), dstar = c(1,1,1))
# hessian(beta =c(1,1,1), y = dat2$VL_1, X = cbind(1, dat2$A2_0, dat2$VL_0), dstar = c(1,1,1))
# set.seed(1)
# bp = beta_prop_mMALA(beta_old = c(1,1,1), epsilon = 0.1, y = dat2$VL_1, X = cbind(1, dat2$A2_0, dat2$VL_0), dstar = c(1,1,1) )
# acc_mMALA(beta_old = c(1,1,1), beta_prop = bp, y = dat2$VL_1, X = cbind(1, dat2$A2_0, dat2$VL_0), epsilon = 0.1, dstar = c(1,1,1))


m0 = glm(VL_1 ~ VL_0 + A2_0, data = dat2, family = binomial())
beta.init = as.numeric(coefficients(m0))
nIter = 2000
epsilon = 0.01
theta.init = 0.5
pstar = c(.5, .5, .5)
nu_0 = 0.05
nu_1 = 5
set.seed(1)

p = 3
params = array(.1, dim = c(nIter, 2 * p + 2))
params[1,5:7] = coefficients(glm(VL_1 ~ VL_0 + A2_0, family = binomial, data = dat2))

########### sample parameters ############
for(it in 2:nIter) {
  # sample theta
  
  ############## sample CD4's #####################
  new.params = sample.new.B(theta.old = params[it-1, 1],
                            gamma.old = params[it-1, 2:(p+1)],
                            beta.old = params[it-1, (p+2):(2*p+1)],
                            epsilon.old = params[it-1, 2*p+2],
                            design = cbind(1, dat2$VL_0, dat2$A2_0),
                            resp = dat2$VL_1 )
  
  params[it, 1] = new.params$theta.new
  params[it, 2:(p+1)] = new.params$gamma.new
  params[it, (p+2):(2*p+1)] = new.params$beta.new
  params[it, 2*p+2] = new.params$epsilon.new
}

plot(params[,6])
acf(params[,5])
