#################################
## Replication codes for simulation study in Bersson and Hoff 2023b
#################################
library(tidyverse)
source("./pred_functions.R")
source("./helper_functions.R")
set.seed(1)
#################################

#################################
## Options
#################################
epsilon = .15
Ks = c(3:150) # num of categories #seq(from =3, to = 150, by = 10) #
Ns = c(10,100,1000)
N.titles = c("N=10","N=100","N=1000")
prob.class = "skewed2" # skewed1,uniform

S  = 500 #5000 # num of data sets in simulation
#################################

#################################
## Run simulation
#################################

out.cov = array(NA,dim = c(length(Ks),length(Ns),
                           4)) # length of output
out.len = array(NA,dim = c(length(Ks),length(Ns),
                           4,# length of output
                           S)
                ) 
dimnames(out.cov)[[1]]=dimnames(out.len)[[1]]=paste0("K,",Ks)
dimnames(out.cov)[[2]]=dimnames(out.len)[[2]]=paste0("N_mult,",round(Ns,1))
dimnames(out.cov)[[3]]=c("oracle","freq","fab_good","fab_bad")
dimnames(out.len)[[3]] = c("fab_good","fab_bad","freq","oracle")
dimnames(out.len)[[4]] = paste0("S",1:S)

alpha.prior = list()



## run for all of these options
for (K.ind in 1:length(Ks)){
  for (N.ind in 1:length(Ns)){
    
  ## get parameters
  K = Ks[K.ind]
  N = Ns[N.ind] 
      
  if (prob.class=="uniform"){
    prob = rep(1,K) * K * 10 
    prob = rdirichlet(prob)
  } else if (prob.class=="skewed1"){
    ## nearly 1/(k/4)
    prob = rep(0,K)
    prob[1:round(K/4)] = 1/round(K/4) * K * 10
    prob = rdirichlet(prob)
    # ## exactly 1/(k/4)
    # prob = rep(0,K)
    # prob[1:round(K/5)] = 1/round(K/5)
    
  } else if (prob.class == "skewed2"){
    prob = rep(1,K)
    # 1/4 big
    prob[1:ceiling(K/4)] = ceiling(K/4):1 * K * 100
    prob = rdirichlet(prob)
  }
  
  oracle.ordering.vec = get.oracle.order.vec(prob)
  
  prec = 100
  # alpha.prior.good = prob*prec
  # if ( prob.class == "uniform"){
  #   alpha.prior.bad = prob*prec
  #   alpha.prior.bad[(round(K/3)+1):K] = 0
  # } else {
  #   alpha.prior.bad = (1-prob)*prec
  # }
  alpha.prior.good = prob*10
  alpha.prior.bad = prob*2
  
      
  ## run code
  ## sim S times
  freq.cov = fab.good.cov = fab.bad.cov = oracle.cov = 0
  freq.len = fab.good.len = fab.bad.len = c()
      
  oracle = oracle_pred(prob,epsilon)
  oracle.len = length(oracle$pred)
  for (s in 1:S){
    
    # obtain data set
    n = rmultinom(1,N,prob) %>% c()
    
    
    # fab pred set
    fab.good = MN_indir_pred(n,alpha.prior.good,epsilon)
    fab.bad = oracle_ordering_pred(n,oracle.ordering.vec,epsilon) # FAB_cat_pred(n,alpha.prior.bad,epsilon) 
    # direct- no prior
    freq = MN_indir_pred(n,rep(0,K),epsilon)
    
    # get cardinality
    fab.good.len = c(fab.good.len, length(fab.good$pred) )
    fab.bad.len = c(fab.bad.len, length(fab.bad$pred) )
    freq.len = c(freq.len, length(freq$pred) )
    
    
    ## is pred covered?
    # predict value from distn'
    y.star = which(rmultinom(1,1,prob)==1)
    
    oracle.cov = oracle.cov + ifelse(y.star%in%oracle$pred,1,0)
    freq.cov = freq.cov + ifelse(y.star %in% freq$pred,1,0)
    fab.good.cov = fab.good.cov + ifelse(y.star %in% fab.good$pred,1,0)
    fab.bad.cov = fab.bad.cov + ifelse(y.star %in% fab.bad$pred,1,0)
    
  }
  
  out.cov[K.ind,N.ind,] = c(oracle.cov,freq.cov,fab.good.cov,fab.bad.cov)/S
  out.len[K.ind,N.ind,,] = rbind(fab.good.len,fab.bad.len,freq.len,length(oracle.len))
      
  }
}
#################################

#################################
## Make plots
#################################
for (N.ind in 1:length(Ns)){
  out.len.Nselec = out.len[,N.ind,,] ## select N = Ns[N.ind]
  plot(NA,NA,
       ylim = range(apply(out.len.Nselec[,1,]/out.len.Nselec[,3,],1,function(j)quantile(j,c(.025,.975))),
                    apply(out.len.Nselec[,2,]/out.len.Nselec[,3,],1,function(j)quantile(j,c(.025,.975)))),
       xlim = range(Ks),
       xlab = "K",
       ylab = "Expected Cardinality Ratio",
       main = paste0(N.titles[N.ind]),
       cex.main = 2,
       cex.axis = 1.5,
       cex.lab = 1.5)
  abline(h=1, lty="dashed")
  for ( k in 2:length(Ks)){
    
    a = out.len.Nselec[k,1,]/out.len.Nselec[k,3,]
    b = out.len.Nselec[k,2,]/out.len.Nselec[k,3,]
    
    points(Ks[k],mean(b),pch = "-",col = "blue",cex=1.5)
    draw_line(Ks[k],mean(b) + c(1,-1)*sd(b)/sqrt(S),alpha("blue",.5))
    points(Ks[k],mean(a),pch = "-", col = "red",cex=1.5)
    draw_line(Ks[k],mean(a) + c(1,-1)*sd(a)/sqrt(S),alpha("red",.5))
  }
  if ( (N.ind == 1) ){
    legend(x = 90, y = 1.3,
           pch="-",
           col=c("red","blue"),
           legend=c(expression(gamma ~ "= 10" ~ theta),expression(o[oracle])),
           cex = 1.5)
  }
}
#################################
