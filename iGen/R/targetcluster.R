targetcluster <- function(u, nclust=41, epsilon=0.04)
{
  set.seed(20)
  
  # define vector for storing clustering results 
  u.ss <- numeric(nclust-1)
  names(u.ss) <- c(2:nclust)


  n <- ncol(u) # of compounds
  m <- nrow(u) # of targets
  
  ## define matrix to store clustring labels
  u.label <- matrix(NA, nrow=m, ncol=nclust-1)
  
  for(ii in 2:41)
  {
    print(ii)
    
    u.kclu <- kmeans(u, ii, nstart=20)
    
    ## assign within cluster distance
    u.ss[ii-1] <- u.kclu$tot.withinss
    ##u.ss[ii-1] <- u.kclu$tot.withinss/ii
    
    ## assign clustering results/labels
    u.label[,ii-1] <- u.kclu$cluster
    
  }
  
  
  ## initialize number of clusters
  nk<- 2
  
  ## choose number of clusters
  while(abs(1-u.ss[nk]/u.ss[nk-1]) >  epsilon){
    nk <- nk+1
  }
  
  
  return(list(label=u.label[,(nk-1)], withinss=u.ss[nk-1], K=nk))
}