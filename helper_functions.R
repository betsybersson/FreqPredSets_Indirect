row_standardize = function(W){
  diag(1/rowSums(W)) %*% W
}
sq_exp_distance_matrix = function(lat,lon,spat.names = 1:length(lat)){
  J = length(lat)
  W = matrix(data=0,nrow=J,ncol=J)
  for ( j in 1:J ){
    
    lon.j = lon[j]
    lat.j = lat[j]
    
    dist.from.j = sapply(1:J,function(k)sqrt((lon.j-lon[k])^2+
                                               (lat.j-lat[k])^2) )
    stand.dist = exp(-dist.from.j^2)
    
    W[j,] = stand.dist
    
  }
  diag(W) = 0
  colnames(W) = rownames(W) = spat.names
  
  return(W)
}

get_nn = function(lat,lon,
                  n.nn = 5,
                  spat.names = 1:length(lat)){
  
  n.counties = length(lat)
  
  W = sq_exp_distance_matrix(lat,lon,
                             spat.names)
  W.stand = row_standardize(W)
  
  Adj = matrix(0,ncol = n.counties, nrow = n.counties)
  for (j.ind in 1:n.counties){
    j = W.stand[j.ind,]
    Adj[j.ind,which(j%in%sort(j,decreasing = T)[1:n.nn])] = 1
  }
  
  Adj.ind = apply(W.stand,1,function(j)which(j%in%sort(j,decreasing = T)[1:n.nn])) %>% t()
  
  return(Adj.ind)
}