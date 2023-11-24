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

michael_cat_pred = function(X,alpha,epsilon,
                            cat_names = 1:length(X)){
  ## michael's form- i think is same thing
  
  K = length(X)
  N = sum(X)
  
  ## test 
  pz.out = array(NA,dim=K)
  names(pz.out) = cat_names
  for ( i in 1:K ){
    
    # Test Z* = i, so add that to our data vector
    Zi = X
    Zi[i] = Zi[i] + 1
    
    pz.out[i] = ((Zi[i] + alpha[i]) >= (Zi + alpha)) %*% Zi / (N+1)
    
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


FAB_cat_pred_standard_CP = function(X,alpha,epsilon,
                             cat_names = 1:length(X)){
  ## X is vector of length K that contains counts from each of the K categories
  ## alpha is vector of length K with dirichlet prior values
  ## epsilon is prediction error rate
  
  K = length(X)
  N = sum(X)
  
  # run conformal algorithm
  # cycle through K categories
  pz.out = array(NA,dim=K)
  names(pz.out) = cat_names
  for ( k in 1:K) {
    X.star = lambda(k,K) # check if this value should be included in prediction set
    conf.scores = MN_ppd(X.star,X,k,alpha) # conf score of X.star
    for ( j in 1:K ){
      X.temp = lambda(j,K) # get conformity score corresponding to this value
      c.temp = MN_ppd(X.temp,X + X.star - X.temp,j,alpha)
      conf.scores = c(conf.scores,rep(c.temp,X[j]))
    }
    
    pz.out[k] = sum(conf.scores <= conf.scores[1])/(N+1)
  }

  epsilon.set = names(which(pz.out>epsilon))
  
  return(list("test_stats" = pz.out,"pred" = epsilon.set))
  
}
FAB_cat_pred_standard_CP_finalcheck = function(X,alpha,epsilon,
                                               cat_names = 1:length(X)){
  ## break X into vector of lambdas
  
  K = length(X)
  N = sum(X)
  
  Y = c()
  for ( k in 1:K){
    Y.temp = matrix(0,nrow = X[k],ncol = K)
    Y.temp[,k] = 1
    Y = rbind(Y,Y.temp)
  }
  Y = rbind(Y,rep(1,K))
  
  
  # run conformal algorithm
  # cycle through N elements for each K
  conf.scores = array(NA,dim=c(N+1))
  pz.out = array(NA,dim=K)
  names(pz.out) = cat_names
  for ( k in 1:K){
    Y[N+1,] = lambda(k,K)
    for ( n in 1:(N+1) ){
      j.ind = which(Y[n,]==1)
      conf.scores[n] = MN_ppd(Y[n,],colSums(Y) - Y[n,],j.ind,alpha)
    }
    pz.out[k] = sum(conf.scores <= conf.scores[N+1])/(N+1)
    
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

#### romano and candes stuff
L.candes = function(pi,tau){
  pi.sort = sort(pi,T)
  pi.sort.cumsum = cumsum(pi.sort)
  pi.sort.cumsum[length(pi.sort.cumsum)] = 1 # numerical thing
  return(min(which(pi.sort.cumsum>=tau)))
}
V.candes = function(pi,tau){
  L.val = L.candes(pi,tau)
  pi.sort = sort(pi,T)
  return(sum(pi.sort[c(1:L.val)]-tau)/pi.sort[L.val])
}
S.candes = function(pi,tau,u){
  ## Eqn 5/6
  ## error rate epsilon
  
  L.val = L.candes(pi,tau)
  pi.sort = sort(pi,T)
  
  set.out=c()
  if (u<=V.candes(pi,tau)){  # with randomization, drop the final group option
    if (length(L.val) == 1){ # set is null set if there's only one option
      set.out = 0
    } else {
      set.out = names(pi.sort)[1:c(L.val-1)]
    }
  } else { # else, keep it
    set.out = names(pi.sort)[1:L.val]
  }
  
  return(set.out)
  
}
E.candes = function(y,pi,u){
  
  check = 0
  tau.out = 0
  tau = 1
  
  # taus.to.check = rev(c(0,cumsum(sort(pi,T)))) ## only have to check at pi.hat change points
  
  while ( check == 0){
    
    S.temp = S.candes(pi,tau,u)
      if(y%in%S.temp){
        tau.out = tau
        tau = tau - .0001
      } else {
        check = 1
      }
    
    # for ( tau in taus.to.check){
    #   S.temp = S.candes(pi,tau,u)
    #   if(y%in%S.temp){
    #     tau.out = tau
    #   } else {
    #     check = 1
    #   }
    # }
        
  }
  
  return(tau.out)
  
}
candes_pred_set = function(n,epsilon){
  
  K = length(n)
  N = sum(n)
  names(n) = c(1:K)
  
  
  y = rep(1:K,times=n)
  # split into 2 subsets
  ind1 = sample(1:N,N/2)
  y1 = y[ind1]
  y2 = y[-ind1]
  

  # sample N+1 us
  U = runif(N+1)
  
  # get pi.hat from data1
  y1.table.temp = table(y1)
  y1.table = rep(0,K)
  y1.table[as.numeric(names(y1.table.temp))] = c(y1.table.temp)
  pi.hat = c(y1.table/(N/2))
  names(pi.hat) = c(1:K)
  
  # compute quantile conformity score on data2 given the trained pi.hat
  E.out = array(NA,c(N/2))
  for ( i in 1:(N/2)){
    E.out[i] = E.candes(y2[i],pi.hat,U[i])
  }
  Q.hat = sort(E.out)[ceiling((1-epsilon)*(1+N/2))]
  
  pred.set = S.candes(pi.hat,Q.hat,U[N+1])
  
  return(pred.set)
  
}
EB_values = function(n.i,X = rep(1,length(n.i)),W,
                     W.stand = row_standardize(W)){
  ## n.i is a vector of length J of observed counts for variable i from all J sites
  ## X covariates
  ## W spatial distance matrix - not row standardized!
  library(spatialreg)
  library(parallel)
  library(spdep)
  
  y = log(n.i+1)
  
  J = length(n.i) # number of sites
  
  df = data.frame("Y" = y, "X" = X)
  
  out = mclapply(1:J,function(j) EB_values_j(j,df,W,W.stand)) %>% 
    unlist()
  
  return(out)
}
EB_values_j = function(j,df,W,W.stand){
  J = nrow(df)
  # get sar mle from data outside of j
  sar.out =  errorsarlm(Y~.-1,data = df[-j,],
                        listw = mat2listw(row_standardize(W[-j,-j])))
  # extract coefs from output
  mu.eb = sar.out$coefficients
  rho.eb = sar.out$lambda
  tau2.eb = sar.out$s2
  # grab spatial covariance
  V = tau2.eb * solve((eye(J) - rho.eb * W.stand) %*% t(eye(J) - rho.eb * W.stand))
  # compute conditional mean given all param
  out = c(df[j,-1] %*% mu.eb + V[j,-j] %*% solve(V[-j,-j]) %*% (df$Y[-j]-mu.eb))
  # exponentiate output to remove log
  out = exp(out) - 1
  ## if less than 0, snap to 0
  out = ifelse(out<0,0,out)
  return(out)
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
conservative_cat_pred = function(X,epsilon,
                               cat_names = 1:length(X)){
  # highest mass region of MHGEOM distribution of X|X+Y
  
  K = length(X)
  N = sum(X)
  
  # get dens for every tested value of predictand Y
  densMHG = sapply(1:K,function(j)dmvhgeom(X,X + lambda(j,K),N+1))
  
  out = HMR(densMHG,epsilon,cat_names)
  
  return(out)
  
}
multinoulli_conservative_cat_pred = function(X,epsilon,
                                 cat_names = 1:length(X)){
  # highest mass region of MHGEOM distribution of X|X+Y
  
  K = length(X)
  N = sum(X)
  
  # get dens for every tested value of predictand Y
  densMHG = sapply(1:K,function(j)(X[j] + 1)/(N+1))
  
  out = HMR(densMHG,epsilon,cat_names)
  
  return(out)
  
}
HMR = function(PROB,epsilon,
               cat.names = 1:length(true.prob)){
  # get highest mass region (HMR) based on probability vector PROB 
  
  K = length(PROB)
  names(PROB) = cat.names
  
  pz.out  = array(NA,dim=K)
  names(pz.out) = cat.names
  for ( i in 1:K ){
    
    pz.out[i] = (PROB[i] >= PROB) %*% PROB
    ## can do cumsum(sort(PROB)) knowing that ordering will change
    
  }
  
  epsilon.set = names(which(pz.out>epsilon))
  
  return(list("test_stats" = pz.out,"pred" = epsilon.set))
  
}
