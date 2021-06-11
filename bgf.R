source('drivers.R')


dat = read.csv("combined2.csv")[,-1]

dat$VL_0 = ifelse(dat$VLOAD_0 < exp(10), 0, 1)
dat$VL_1 = ifelse(dat$VLOAD_1 < exp(10), 0, 1)

#attach(dat.f1)
set.seed(2)
nIter = 1000

nu_1 = 5
nu_0 = 0.3

p = 2 # covariate model design matrix dimension
pY = 4 # outcome model design matrix dimension

# array of parameters for each intermediate variable
# theta (1)/ gamma (p)/ beta (p)/ sigsq (1)

nCov = 2 # number of time-indexed covariates
params = array(.1, dim = c(2, nIter, 2 * p + 2)) # intermediate variables
paramsB = array(.1, dim = c(1, nIter, 2 * p + 2)) # params for binary variable (vload)
paramsY = array(.1, dim = c(nIter, 2 * pY + 2))

########### sample parameters ############
for(it in 2:nIter) {
  # sample theta
  
  ############## sample VLOAD's #####################
  new.params = sample.new.B(theta.old = paramsB[1, it-1, 1],
                            gamma.old = paramsB[1, it-1, 2:(p+1)],
                            beta.old = paramsB[1, it-1, (p+2):(2*p+1)],
                            epsilon.old = params[1, it-1, 2*p+2],
                            design = cbind( dat$VL_0, dat$A2_0),
                            resp = dat$VL_1 )
  
  paramsB[1, it, 1] = new.params$theta.new
  paramsB[1, it, 2:(p+1)] = new.params$gamma.new
  paramsB[1, it, (p+2):(2*p+1)] = new.params$beta.new
  paramsB[1, it, 2*p+2] = new.params$epsilon.new
  
  
  ############## sample CD4's #####################
  new.params = sample.new.N(theta.old = params[1, it-1, 1],
                          gamma.old = params[1, it-1, 2:(p+1)],
                          beta.old = params[1, it-1, (p+2):(2*p+1)],
                          sigsq.old = params[1, it-1, 2*p+2],
                          design = cbind(log(dat$CD4N_0), dat$A2_0),
                          resp = log(dat$CD4N_1) )
  
  params[1, it, 1] = new.params$theta.new
  params[1, it, 2:(p+1)] = new.params$gamma.new
  params[1, it, (p+2):(2*p+1)] = new.params$beta.new
  params[1, it, 2*p+2] = new.params$sigsq.new
  
  ############## sample CD8's #####################
  new.params = sample.new.N(theta.old = params[2, it-1, 1],
                          gamma.old = params[2, it-1, 2:(p+1)],
                          beta.old = params[2, it-1, (p+2):(2*p+1)],
                          sigsq.old = params[2, it-1, 2*p+2],
                          design = cbind(log(dat$CD8N_0), dat$A2_0),
                          resp = log(dat$CD8N_1) )
  
  params[2, it, 1] = new.params$theta.new
  params[2, it, 2:(p+1)] = new.params$gamma.new
  params[2, it, (p+2):(2*p+1)] = new.params$beta.new
  params[2, it, 2*p+2] = new.params$sigsq.new
  
  ############## sample Y's #######################
  pY = 4
  new.params = sample.new.N(theta.old = paramsY[it-1, 1],
                          gamma.old = paramsY[it-1, 2:(pY+1)],
                          beta.old = paramsY[it-1, (pY+2):(2*pY+1)],
                          sigsq.old = paramsY[it-1, 2*pY+2],
                          design = cbind(log(dat$CD4N_1), log(dat$CD8N_1), dat$A2_0, dat$A2_1),
                          resp = log(dat$CD4N_2) )
  
  paramsY[it, 1] = new.params$theta.new
  paramsY[it, 2:(pY+1)] = new.params$gamma.new
  paramsY[it, (pY+2):(2*pY+1)] = new.params$beta.new
  paramsY[it, 2*pY+2] = new.params$sigsq.new
  
}


#### GF using posterior predictive density

# completed dataset for each set of parameters
# RANDOM G
# augmented data set
n = nrow(dat)

dat$y11 = NA
dat$y00 = NA
dat$t1 = ifelse(dat$A2_0 == 1 & dat$A2_0 == 1, 1, 0)
dat$t0 = ifelse(dat$A2_0 == 0 & dat$A2_0 == 0, 1, 0)

# first, let's do this with the correct G
row = 100
ndraws = 900

diffs = array(0, dim = ndraws)

for(draw in 1:ndraws) {
  
  if ( params[1, row + draw, 3] == 1 ) { # if there is an edge from treatment to covariate, impute
    
    mx1 = cbind(log(dat$CD4N_0), 1) %*% params[1, row + draw, c(4,5)]  # XB (treated at time 1)
    X2draw1 = rnorm(n, mean = mx1, sd = params[1, row + draw, 6]) # fixed treatment to 1
    X2draw1 = ifelse(dat$A2_0 == 1, log(dat$CD4N_0), X2draw1)
    
    
    mx0 = cbind(log(dat$CD4N_0),0) %*% params[1, row + draw, c(4,5)]  # XB (treated at time 1)
    X2draw0 = rnorm(n, mean = mx0, sd = params[1, row + draw, 6]) # fixed treatment to 1
    X2draw0 = ifelse(dat$A2_0 == 0, log(dat$CD4N_0), X2draw0)
    
    
  } else { 
    
    X2draw1 = log(dat$CD4N_0)
    X2draw0 = log(dat$CD4N_0)
    
  }
  
  if ( params[2, row + draw, 3] == 1 ) { # if there is an edge from treatment to covariate, impute
    
    mz1 = cbind(log(dat$CD8N_0), 1) %*% params[2, row + draw, c(4,5)]  # XB (treated at time 1)
    Z2draw1 = rnorm(n, mean = mz1, sd = params[2, row + draw, 6]) # fixed treatment to 1
    Z2draw1 = ifelse(dat$A2_0 == 1, log(dat$CD8N_0), Z2draw1)
    
    
    mz0 = cbind(log(dat$CD8N_0), 0) %*% params[1, row + draw, c(4,5)]  # XB (treated at time 1)
    Z2draw0 = rnorm(n, mean = mz0, sd = params[1, row + draw, 6]) # fixed treatment to 1
    Z2draw0 = ifelse(dat$A2_0 == 0, log(dat$CD8N_0), Z2draw0)
    
    
  } else { 
    
    Z2draw1 = log(dat$CD8N_0)
    Z2draw0 = log(dat$CD4N_0)
    
  }
  
  my1 = cbind(X2draw1, Z2draw1, 1, 1) %*% paramsY[(row+draw),c(6:9)] # always treated
  Ydraw1 = rnorm(n, mean = my1, sd = paramsY[(row+draw), 10])
  
  my0 = cbind(X2draw0, Z2draw0, 0, 0) %*% paramsY[(row+draw),c(6:9)] # always treated
  Ydraw0 = rnorm(n, mean = my0, sd = paramsY[(row+draw), 10])
  
  
  dat$y11 = ifelse(dat$t1 == 1, log(dat$CD4N_2), Ydraw1)
  dat$y00 = ifelse(dat$t0 == 1, log(dat$CD4N_2), Ydraw0)
  
  #print(sum(dat$y11))
  diffs[draw] = sum(dat$y11 - dat$y00)/n
}



