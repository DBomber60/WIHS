library(mvtnorm)

# log probability of beta given data, gamma. elements of dstar are diagonal elements of D^(-1)
logpi = function(beta,y,X, dstar){
  eta = X%*%beta;
  return(sum(y*eta - log(1+exp(eta))) - .5 * sum( dstar * beta^2 ) )
}

# compute gradient based on a N(0,D) prior on beta
grad = function(beta, y, X, dstar) {
  eta = X%*%beta;
  return(-dstar * beta + t(X) %*% (y - plogis(eta)))
}

# compute inverse hessian based on a N(0,D) prior on beta
# in calderman parlance, this is G^-1
hessian = function(beta,y,X, dstar) {
  eta = X%*%beta;
  phat = plogis(eta)
  W = as.numeric(phat* (1-phat)) # Var(y)
  C = solve(diag(as.numeric(dstar)) + t(as.numeric(phat* (1-phat)) * X ) %*% X)
  return(C)
}

# mMALA beta proposal
beta_prop_mMALA = function(beta_old, epsilon, y, X, dstar) {
  p = dim(X)[2]
  C = hessian(beta_old, y, X, dstar)
  rootC = chol(C)
  return( beta_old + .5 * epsilon^2 * C %*% grad(beta_old, y, X, dstar) + epsilon * rootC %*% rnorm(p) )
}

# MMALA acceptance probability
acc_mMALA = function(beta_old, beta_prop, y, X, epsilon, dstar) {
  a = exp(logpi(beta_prop,y,X, dstar) - logpi(beta_old,y,X, dstar) )
  
  Ginv_old = hessian(beta_old, y, X, dstar)
  grad_old = grad(beta_old, y, X, dstar)
  
  Ginv_new = hessian(beta_prop, y, X, dstar)
  grad_new = grad(beta_prop, y, X, dstar)
  
  new_given_old = dmvnorm(as.numeric(beta_prop), mean= beta_old + .5 * epsilon^2 * Ginv_old %*% grad_old,
                          sigma = epsilon^2 * Ginv_old, log = T)
  
  old_given_new = dmvnorm(as.numeric(beta_old), mean= beta_prop + .5 * epsilon^2 * Ginv_new %*% grad_new,
                          sigma = epsilon^2 * Ginv_new, log = T)
  
  b = exp(old_given_new - new_given_old)
  return(a * b)
}





# sample new set of parameter values for a normally distributed response
sample.new.B = function(theta.old, gamma.old, beta.old, epsilon.old, design, resp) {
  
  p = ncol(design)
  # sample new theta
  # assumption: U(0,1) prior on theta/ binomial distribution on gamma (beta - binomial)
  
  theta.new = rbeta(1, shape1 = sum(gamma.old) + 1, shape2 = p - sum(gamma.old) + 1)
  
  # sample new gamma 
  # assumption: independent elements, fixed nu_0, nu_1
  
  p1 = dnorm(beta.old, sd = sqrt(nu_1)) * theta.new
  p0 = dnorm(beta.old, sd = sqrt(nu_0)) * (1-theta.new)
  gamma.new = rbinom(p, 1, prob = p1/(p1 + p0))  

  # sample new beta
  dstar = ifelse(gamma.new == 1, 1/nu_1, 1/nu_0)
  beta.prop = beta_prop_mMALA(beta.old, epsilon.old, resp, design, dstar)
  #print(beta.prop)
  
  # accept it?
  al = acc_mMALA(beta_old = beta.old, beta_prop = beta.prop, y = resp, X = design, epsilon = epsilon.old, dstar = dstar)
  print(al)
  Acc = min(1,al)
  if (runif(1)<=Acc){
    beta.new = beta.prop;
  } else {beta.new = beta.old}
  
  epsilon.new = exp( log(epsilon.old) + (1/it^0.7)*(Acc - 0.4) )
  print(epsilon.new)
  
  # assume scale mixture prior and normal likelihood for response with sigma = sigma^2 * I
  #M = (t(design) %*% design)*sigsq.old^-1 + diag(ifelse(gamma.new == 1, 1/nu_1, 1/nu_0)) 
  #Minv = solve(M)
  #meanvec = sigsq.old^-1 * Minv %*% t(design) %*% resp
  #beta.new = rmvnorm(1, mean = meanvec, sigma = Minv)
  
  # sample new sigsq
  # assumptions: IG(.5, .5) prior on sigsq
  #ss = sum ( (resp - design %*% array(beta.new, dim=p) )^2 )
  #n = length(resp)
  #sigsq.new = 1/rgamma(1, (n+1)/2, (ss+1)/2)
  
  
  return(list(theta.new=theta.new, 
              gamma.new=gamma.new, 
              beta.new=beta.new,
              epsilon.new=epsilon.new))
}

