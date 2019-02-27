

#registerDoParallel(5)  
gc()
library("rlang")
library('igraph')
library("bayess")
library('rlist')
library('MASS')
library('matrixcalc')
library('MCMCpack')
library('inline')
library('RcppArmadillo')
library('Rcpp')
#library("rstan")
library("Boom")
#library("profvis")
library("mclust")
library("plot3D")
library("geoR")
devtools::find_rtools()

source('MAIN_FUNCTION.R')
result_all<-list()

for(seed in 1:10){
 
  set.seed(seed)
  
  
  
  
  c=40
  g_full <- graph.lattice(length=c,dim=2)
  n=c^2
  #V(g_full)$index=1:n
  net=get.adjacency(g_full,attr=NULL)
  K=5
  rho=5
  label_T0<-pottshm(ncol=K,niter=1,c,m=c,beta=rho)
  label_T1<-pottshm(ncol=K,niter=1,c,m=c,beta=rho)
  cc=c/4
  
  label_T0[1:cc,]<-1
  label_T0[(cc+1):(cc*2),]<-2
  label_T0[(cc*2+1):(cc*3),]<-3
  label_T0[(cc*3+1):c,]<-4
  Matrix_list=lapply(1:K, function(k) rwish(30, (k+1)*diag(3)))
  Matrix_list[[K]]<-rwish(10, 1.5*diag(3))
  label_T1[1:cc,]<-1
  label_T1[(cc+1):(cc*2),]<-2
  label_T1[(cc*2+1):(cc*3),]<-3
  label_T1[(cc*3+1):c,]<-4
  
  label_T1[(cc+1):(cc*2),(cc+1):(cc*2)]<-5
  
  
  Data=list()
  N_T0=5
  N_T1=5
  T_index=c(rep(1,N_T0),rep(2,N_T1))
  
  Data_T0<-lapply(1:N_T0, function(i) lapply(1:c^2,function(s) riwish(3,Matrix_list[[label_T0[s]]])))
  Data_T1<-lapply(1:N_T1, function(i) lapply(1:c^2,function(s) riwish(3,Matrix_list[[label_T1[s]]])))
  Data<-append(Data_T0,Data_T1)
  par(mfrow=c(1,2))
  image(label_T0)
  image(label_T1)
  
  iter=5000
  burnin=3000
  p=3
  n=c^2
  K=10
  N=N_T0+N_T1
  no.T=2
  Ca=3
  Cb=3
  adj_list<- apply(net==1,1,which)
  S=diag(2)*0.1
  rho_alpha=0.3
  rho_beta=2
  
  
  try(result<-Potts_Bayesian_Semi_Multi(Data,net,n,p=3,N,
                                    K=10,T_index,iter=5000,
                                    LU_alpha=c(2,10), 
                                    LU_beta=c(2,10), 
                                    LU_xi=c(0,1),
                                    rho_alpha=4,
                                    rho_beta=4,
                                    rho_xi=0.5),TRUE)
  result_in<-duplicate(result, shallow = FALSE)
  result_all[[seed]]<-result_in
 
  
  
  
  
  
}

save.image("K10.Rdata")




