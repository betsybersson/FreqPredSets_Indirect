MN_indir_pred = function(X,alpha,epsilon,
                        cat_names = 1:length(X)){
  ## X is vector of length K that contains counts from each of the K categories
  ## alpha is vector of length K with dirichlet prior values
  ## epsilon is prediction error rate
  
  K = length(X)
  if ( K != length(alpha)){
    ## error!
  }
  N = sum(X)
  
  
  # for each category, obtain xl + alphal
  c_i = X + alpha
  
  # for each category, obtain c_nPLUS1 in that case
  c_nPLUS1 = X + 1 + alpha
  
  ## test 
  pz.out = array(NA,dim=K)
  names(pz.out) = cat_names
  for ( k in 1:K ){
    
    # pz = c(c_nPLUS1[k] >= c_i + lambda(k,K)) %*% X + 1  # equivalent to below
    pz = c(c_nPLUS1[k] >= c_i) %*% X + 1 
    # normalize
    pz = pz/(N+1)
    
    pz.out[k] = pz
    
  }
  
  epsilon.set = names(which(pz.out>epsilon))
  
  return(list("test_stats" = pz.out,"pred" = epsilon.set))
  
  
}
get.oracle.order.vec = function(PROB,cat_names = 1:length(PROB)){
  ### no ties!!!
  order.prob = order(PROB)
  o.out = sapply(PROB,function(j)which(j == PROB[order.prob]))
  names(o.out) = cat_names
  return(o.out)
}
oracle_ordering_pred = function(X,o,epsilon,
                        cat_names = 1:length(X)){
  ## X is vector of length K that contains counts from each of the K categories
  ## epsilon is prediction error rate
  
  K = length(X)
  if ( K != length(o)){
    ## error!
  }
  N = sum(X)
  
  ## test 
  pz.out = array(NA,dim=K)
  names(pz.out) = cat_names
  for ( k in 1:K ){
    
    # pz = c(c_nPLUS1[k] >= c_i + lambda(k,K)) %*% X + 1  # equivalent to below
    pz = c(o[k] >= o) %*% X + 1 
    # normalize
    pz = pz/(N+1)
    
    pz.out[k] = pz
    
  }
  
  epsilon.set = names(which(pz.out>epsilon))
  
  return(list("test_stats" = pz.out,"pred" = epsilon.set))
  
  
}

bayesian_cat_pred = function(X,alpha,epsilon,
                             cat_names = 1:length(X)){
  # bayesian posterior predictive set

  K = length(X)
  N = sum(X)
  
  ## test 
  pz.out  = array(NA,dim=K)
  names(pz.out) = cat_names
  for ( i in 1:K ){
    
    pz.out[i] = ((X[i] + alpha[i]) >= (X + alpha)) %*% (X + alpha) / (N+sum(alpha))
    
  }
  
  epsilon.set = names(which(pz.out>epsilon))
  
  return(list("test_stats" = pz.out,"pred" = epsilon.set))
  
}
oracle_pred = function(true.prob,epsilon,
                       cat.names = 1:length(true.prob)){
  # true oracle pred region from probability vector true.prob 
  
  K = length(true.prob)
  names(true.prob) = cat.names
  
  pz.out  = array(NA,dim=K)
  names(pz.out) = cat.names
  for ( i in 1:K ){
    
    pz.out[i] = (true.prob[i] >= true.prob) %*% true.prob
    
  }
  
  epsilon.set = names(which(pz.out>epsilon))
  
  return(list("test_stats" = pz.out,"pred" = epsilon.set))
  
}


lambda = function(k,K){
  out = rep(0,K)
  out[k] = 1
  return(out)
}
MN_ppd = function(y_star,y,k,alpha){
  # MN for K categories; each entry should be vector of dimension K
  # y_star is predictive vector- want to obtain density at this value
  # y is sum of observations of each cateogry
  # alpha is dirichlet prior hyperparameter
  
  alpha_tick = y+alpha
  
  # ppd = lgamma(sum(y_star+1)) - (sum(lgamma(y_star+1))) +
  #   sum(lgamma(y_star + alpha_tick)) + lgamma(sum(alpha_tick)) -
  #   (lgamma(sum(y_star+alpha_tick)) + sum(lgamma(alpha_tick)))
  # ppd = exp(ppd)
  
  # simplifies to this for multinoulli prediction
  ppd = alpha_tick[k]/(sum(alpha) + sum(y))
  
  return(ppd)
  
}

rdirichlet = function(alpha){
  Y = sapply(alpha,function(j)rgamma(1,j,1))
  
  return(Y/sum(Y))
}


#### functions to get dirichlet multinomial mle
##########################################################
## obtain MLE of Polya (dirichlet- multinmonial) distribution
# D: JxK matrix of counts; each row is a sample from a MN distribution with K categories
# init: If NA, use method moment matching procedure to obtain good init values
# method: Choose from "Newton_Raphson", "fixed_point", "separate", "precision_only"
# epsilon: convergence diagnostic
# print_progress: if T, print progress to screen
##########################################################
polyaMLE = function(D, init = NA, method = "Newton_Raphson",
                    epsilon = .0001,
                    print_progress = F) {
  library(parallel)
  
  ## if a column is all 0s, remove it from analysis and set prior to 0
  zeros.ind = which(colSums(D)==0)
  if (length(zeros.ind)>0){
    not.zeros.ind = c(1:ncol(D))[-zeros.ind]
    D.orig = D
    D = D[,-zeros.ind]
  }
  
  if (is.na(init)){
    init = init.fn(D) 
  }
  
  # grab helpers
  Nj = rowSums(D)
  J = nrow(D)
  K = ncol(D)
  
  # some initialization steps
  alpha_t = init
  if (method == "separate"| method == "precision_only"){
    theta0_t = colMeans(row_standardize(D))
    alpha0_t = mean(rowSums(D))
    nu.fn = function(nn,alpha.star){
      sum(alpha.star * (digamma(nn + alpha.star) - digamma(alpha.star)))
    }
  }
  converged = F
  # iterate until convergence
  while (!converged) {
    
    if (method == "Newton_Raphson"){
      # get gradient and inverse Hessian
      g = polyaGradient(D,alpha_t,Nj,K)
      H = polyaHessian(D,alpha_t,Nj,K)$H
      H.inv = solve(H) # should be able to reduce comp time a bit more than this
      
      # update alpha 
      alpha_tP1 = alpha_t - c(H.inv %*% g)
    } else if (method == "fixed_point") {
      # get helper
      alpha0t = sum(alpha_t)
      gprime.k = unlist(mclapply(1:K,function(k) sum(digamma(D[,k]+alpha_t[k]) -
                                                       digamma(alpha_t[k]))/
                                   sum(digamma(Nj + alpha0t) - digamma(alpha0t))
      ))
      
      # update alpha
      alpha_tP1 = alpha_t*gprime.k
    } else if (method == "separate"){
      # update precision- alpha0
      gprime.k = sum(unlist(mclapply(1:K,function(k) 
        sum(theta0_t[k] * digamma(D[,k]+ alpha0_t * theta0_t[k]) -
              theta0_t[k] * digamma(alpha0_t * theta0_t[k])))))/
        sum(digamma(Nj + alpha0_t) - digamma(alpha0_t))
      # update alpha0
      alpha0_tP1 = alpha0_t*gprime.k
      
      ## update theta0
      theta0k =  unlist(mclapply(1:K,function(k)nu.fn(D[,k],alpha0_tP1*theta0_t[k])))
      theta0_tP1 = theta0k/sum(theta0k)
      
      # compute alphat for convergence check and output
      alpha_tP1 = theta0_tP1*alpha0_tP1
      
      # save updated values
      alpha0_t = alpha0_tP1
      theta_t = theta0_tP1
      
    }else if (method == "precision_only"){
      # update precision- alpha0
      gprime.k = sum(unlist(mclapply(1:K,function(k) 
        sum(theta0_t[k] * digamma(D[,k]+ alpha0_t * theta0_t[k]) -
              theta0_t[k] * digamma(alpha0_t * theta0_t[k])))))/
        sum(digamma(Nj + alpha0_t) - digamma(alpha0_t))
      # update alpha0
      alpha0_tP1 = alpha0_t*gprime.k
      
      # compute alphat for convergence check and output
      alpha_tP1 = theta0_t*alpha0_tP1
      
      # save updated values
      alpha0_t = alpha0_tP1
    }
    
    
    # check convergence
    converged = ifelse(all(abs(alpha_tP1-alpha_t)<epsilon),T,F)
    alpha_t = ifelse(alpha_tP1<0,1e-10,alpha_tP1) # enforce bounded below by 0 constraint
    
    if (print_progress){
      print(paste0(c("alpha:",round(alpha_t,2)), collapse= " "))
    }
    
  }
  
  # if some species are all 0, put prior of 0 back in
  if (length(zeros.ind)>0){
    alpha.out = array(0,dim=ncol(D.orig))
    alpha.out[not.zeros.ind] = alpha_t
    alpha_t = alpha.out
  }
  
  return(alpha_t)
  
}
##########################################################
## obtain gradient of Polya (dirichlet- multinmonial) distribution;
## output: gradient: first derivative vector of length K
# D: JxK matrix of counts; each row is a sample from a MN distribution with K categories
# alpha: values of prior alpha
# Nj: rowsums of D
# K: dimension of MN (i.e. number of columns in D)
##########################################################
polyaGradient = function(D,alpha,Nj,K){
  ## helpers
  alpha0 = sum(alpha)
  
  g_k = unlist(
    mclapply(1:K,function(k) sum(digamma(alpha0) -
                                   digamma(alpha0 + Nj) +
                                   digamma(D[,k] + alpha[k]) -
                                   digamma(alpha[k])))
  )
  
  return(g_k)
}
##########################################################
## obtain Hessian of Polya (dirichlet- multinmonial) distribution;
## output: Hessian: second derivative matrix of length KxK
# D: JxK matrix of counts; each row is a sample from a MN distribution with K categories
# alpha: values of prior alpha
# Nj: rowsums of D
# K: dimension of MN (i.e. number of columns in D)
##########################################################
polyaHessian= function(D,alpha,Nj,K){
  ## helpers
  alpha0 = sum(alpha)
  
  z = sum(trigamma(alpha0) - 
            trigamma(alpha0 + Nj))
  
  q.diag = unlist(mclapply(1:K,function(k)
    sum(trigamma(D[,k] + alpha[k]) -
          trigamma(alpha[k]))
  ))
  Q = diag(q.diag)
  
  H = Q + z
  
  return(list("q.diag" = q.diag,"z" = z,"H" = H))
}
##########################################################
## obtain inital guess of MLE of Polya (dirichlet- multinmonial) distribution;
## use method of moment matching as in Minka (2000)
# D: JxK matrix of counts; each row is a sample from a MN distribution with K categories
##########################################################
init.fn = function(D){
  p = t(apply(D,1,function(l)l/sum(l)))
  exp.p = colMeans(p)
  exp2.p = colMeans(p*p)
  ok = exp.p>0
  sum.as = (exp.p[ok] - exp2.p[ok]) / (exp2.p[ok] - exp.p[ok]^2)
  init = exp.p * median(sum.as)
  return(init)
}
##########################################################
## multivariate hypergeometric density
## evaluate mulviaraiate hypergeometric density at vector X
# X: a realisation; vector of length m
# K: vector of counts of success states; vector of length m
# N: population size
# n: number of draws, sum of X
##########################################################
dmvhgeom = function(X,K,N,
                    n = sum(X)){
 m = length(X)
 out = prod(sapply(1:m,function(mm)choose(K[mm],X[mm]))) /
   choose(N,n)
 return(out)
}

