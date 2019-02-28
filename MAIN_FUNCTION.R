## Copyright (C) 2018 Zhou Lan (zlan@ncsu.edu) - All Rights Reserved
## You may use, distribute and modify this code under the
## * terms of the  GPL-2. GPL-3, 
## *
## * You should have received a copy of the GPL-2. GPL-3 license with
## * this file. If not, please write to: zlan@ncsu.edu, or visit: https://www4.stat.ncsu.edu/~zlan/Publication.html
## */





library(rlang)
library('igraph')
library("bayess")
library('rlist')
library('MASS')
library('matrixcalc')
library('MCMCpack')
library('inline')
library('RcppArmadillo')
library('Rcpp')
library("rstan")
library("Boom")
#library("profvis")
library("mclust")
library("plot3D")
library("geoR")
#devtools::find_rtools()





update.group.level<-function(nn,kk){
  
  if(kk %in% unique_label){
    result=rho_beta*sum(label_group[adj_list_input[[nn]],t]==kk)+rho_alpha*sum(label_sub[nn,]==kk)
  }else{
    result= rho.alpha*sum(label_sub[nn,]==kk)
  }
  if(is.na(result)){
    result=0
  }
  return(result)
  
}


log.diwish.partial<-function(data,dof,matrix){
  
  0.5*dof*log(det(matrix))-0.5*sum(diag((data)%*%matrix))
  
}

logdet<-function(M){
  return(log(det(M)))
}





CPP_dof_iw.body<-'
RNGScope scope;
IntegerVector n_ll = n_l;
int n=n_ll.size();
NumericVector nu=dof;
NumericVector dimension=p;
NumericVector logdetsumsq=Matrix_list_log_det_input;
NumericVector logdetSigma=logdetData;
NumericVector lmg_V=lmgvalue;
IntegerVector g_subject=ll_sub;

NumericVector mysum(1);

for(int i=0; i<n; i++) {

mysum[0]=mysum[0]-0.5*nu[0]*logdetsumsq[g_subject[i]-1]-0.5*(nu[0]+dimension[0]+1)*logdetSigma[i]-0.5*nu[0]*dimension[0]*log(2)-lmg_V[0];
}


NumericVector mysum_output(clone(mysum));
return(wrap(mysum_output));


'

CPP_dof_iw<-cxxfunction(signature(dof="numeric",p="numeric",Matrix_list_log_det_input="numeric",logdetData="numeric",ll_sub="integer",lmgvalue="numeric",n_l="integer"), body = CPP_dof_iw.body,plugin="RcppArmadillo")



CPP_new.matrix.body<-'
RNGScope scope;
IntegerVector G_subject = ll_s;  
IntegerVector Kll = K_l;  
int K = Kll.size();
int n = G_subject.size();
Rcpp::List inv_Data_clist(inv_Data_list);
Rcpp::List rrr(initial_list);
for (int k=0; k<K; k++){
for (int i=0; i<n; i++){
if (Kll[k]==G_subject[i]){
SEXP data_inv =  inv_Data_clist[i];
arma :: mat data_inv_c=Rcpp :: as < arma :: mat >(data_inv);

SEXP initial =  rrr[k];
arma :: mat op=Rcpp :: as < arma :: mat >(initial);

arma :: mat input=op+data_inv_c;
rrr[k]=input;
}else{


SEXP initial =  rrr[k];
arma :: mat op=Rcpp :: as < arma :: mat >(initial);

arma :: mat input=op;
rrr[k]=input;

}



}
}


Rcpp::List rrr_output(clone(rrr));
return(wrap(rrr_output));

'

CPP_new.matrix<-cxxfunction(signature(ll_s="integer",K_l="integer",inv_Data_list="List",initial_list="List"),body = CPP_new.matrix.body,plugin="RcppArmadillo")



CPP_new.dof.body<-'
RNGScope scope;
IntegerVector G_subject = ll_s;  
NumericVector Kll = K_l;  
int K = Kll.size();
int BN=G_subject.size();

NumericVector dof_new_temp(K);

for (int k=0; k<K; k++){
for (int bn=0; bn<BN; bn++){
if (Kll[k]==G_subject[bn]){
dof_new_temp[k]=dof_new_temp[k]+1;}
}
}


NumericVector dof_new_temp_output(clone(dof_new_temp));
return(wrap(dof_new_temp_output));


'

CPP_new.dof<-cxxfunction(signature(ll_s="integer",K_l="integer"),body = CPP_new.dof.body,plugin="RcppArmadillo")


V2_CPP_rho_alpha.body<-'
RNGScope scope;
IntegerVector n_ll = n_l;
int n=n_ll.size();
NumericVector k_ll = k_l;
int K=k_ll.size();
IntegerVector NT_ll = NT_l;
int NT=NT_ll.size();
IntegerVector g_subject_t=ll_sub_t;


NumericMatrix rr(n,K);
double PHIA = Rcpp::as<double>(phi_alpha);
for(int k=0; k<K; k++){
for(int i=0; i<n; i++) {
rr(i,k)=0;
for(int NNT=1; NNT<(NT+1); NNT++){
if(g_subject_t[i*NT+NNT-1] == k_ll[k]){
rr(i,k)=rr(i,k)+PHIA;}
} 
}
}

NumericVector rr_output(clone(rr));
return(wrap(rr_output));

'

V2_CPP_rho_alpha<-cxxfunction(signature(ll_sub_t="integer",k_l="integer",NT_l="integer",n_l="integer",phi_alpha="integer"), body = V2_CPP_rho_alpha.body,plugin="RcppArmadillo")


CPP_indicator.body<-'
RNGScope scope;
NumericVector ss_rho_c=ss_rho_all_input;
IntegerVector G_g1 = ll_g1;  
IntegerVector G_g2 = ll_g2;  
IntegerVector n_ll = n_l;
Rcpp::List alist(adj_list_input);
int n=n_ll.size();
IntegerVector Kll = K_l;  
int K = Kll.size();
double PHIB = Rcpp::as<double>(phi_beta);


NumericVector rr(n);
NumericVector prior_g1(K);
NumericVector prior_g2(K);
NumericVector prior_sub(K);
NumericVector den_g1(1);
NumericVector den_g2(1);
NumericVector den_sub(1);
NumericVector den(1);
NumericVector no(1);
NumericVector ratio(1);
NumericVector U(1);
NumericVector rtemp(1);


for(int i=0; i<n; i++) {



SEXP ll = alist[i];
Rcpp::IntegerVector nei(ll);
int m=nei.size();  
for(int k=0; k<K; k++){
prior_g1[k]=-100;
prior_g2[k]=-100;
prior_sub(k)=-100;
}


for(int j=0; j<m; j++){
prior_g1[G_g1[nei[j]-1]-1]=prior_g1[G_g1[nei[j]-1]-1]+PHIB;
}

den_g1[0]=0;
for(int k=0;k<K; k++){
den_g1[0]=den_g1[0]+exp(prior_g1[k]);
}




for(int j=0; j<m; j++){
prior_g2[G_g2[nei[j]-1]-1]=prior_g2[G_g2[nei[j]-1]-1]+PHIB;
}
den_g2[0]=0;
for(int k=0;k<K; k++){
den_g2[0]=den_g2[0]+exp(prior_g2[k]);
}



for(int k=0;k<K; k++){
prior_sub(k)=ss_rho_c[k*n-1+i];
}
den_sub[0]=0;
for(int k=0;k<K; k++){
den_sub[0]=den_sub[0]+exp(prior_sub[k]);
}

den[0]=den_g1[0]*den_g2[0]*den_sub[0];


no[0]=0;
for(int k=0;k<K; k++){
no[0]=no[0]+exp(prior_g1[k]+prior_g2[k]+prior_sub[k]);
}


ratio[0]=no[0]/den[0];

U = runif(1,0,1);
rr[i]=1;
if(U[0]<ratio[0]){
rr[i]=0;
}




}




NumericVector rr_output(clone(rr));
return(wrap(rr_output));


'


CPP_indicator<-cxxfunction(signature(ss_rho_all_input="numeric",ll_g1="integer", ll_g2="integer",adj_list_input="List",n_l="integer",K_l="integer",phi_beta="numeric"), body = CPP_indicator.body,plugin="RcppArmadillo")





CPP_update.group.level.body<-'
RNGScope scope;

IntegerVector G_group = ll_g;  
IntegerVector Kll = K_l;  
NumericVector s_sum_rho = ss_rho_input;  
Rcpp::List alist(adj_list_input);
int n = G_group.size();
int K = Kll.size();
double PHIB = Rcpp::as<double>(phi_beta);

NumericVector prior(K);
NumericVector U(1);
NumericVector CDF(K);


for(int i=0; i<n; i++) {



SEXP ll = alist[i];
Rcpp::IntegerVector nei(ll);
int m=nei.size();  
for(int k=0; k<K; k++){
prior[k]=-100;
CDF[k]=0;
}
for(int j=0; j<m; j++){
prior[G_group[nei[j]-1]-1]=prior[G_group[nei[j]-1]-1]+PHIB;
}



CDF[0]=exp(prior[0]+s_sum_rho[i-1]);



for(int k=1; k<K; k++){

CDF[k]=CDF[k-1]+exp(prior[k]+s_sum_rho[k*n-1+i]);

}



U = CDF[K-1]*runif(1,0,1);
G_group[i]=1;
for(int k=1; k<K; k++){
if(U[0]>CDF[k-1]){G_group[i]=k+1;}
}



}


IntegerVector G_group_output(clone(G_group));
return(wrap(G_group_output));


'

CPP_update.group.level<-cxxfunction(signature(ll_g="integer",ss_rho_input="numeric" ,adj_list_input="List",K_l="integer",phi_beta="numeric"), body = CPP_update.group.level.body,plugin="RcppArmadillo")





V2_CPP_rho_alpha.body<-'
RNGScope scope;
IntegerVector n_ll = n_l;
int n=n_ll.size();
NumericVector k_ll = k_l;
int K=k_ll.size();
IntegerVector NT_ll = NT_l;
int NT=NT_ll.size();
IntegerVector g_subject_t=ll_sub_t;


NumericMatrix rr(n,K);
double PHIA = Rcpp::as<double>(phi_alpha);
for(int k=0; k<K; k++){
for(int i=0; i<n; i++) {
rr(i,k)=0;
for(int NNT=1; NNT<(NT+1); NNT++){
if(g_subject_t[i*NT+NNT-1] == k_ll[k]){
rr(i,k)=rr(i,k)+PHIA;}
} 
}
}

NumericVector rr_output(clone(rr));
return(wrap(rr_output));

'

V2_CPP_rho_alpha<-cxxfunction(signature(ll_sub_t="integer",k_l="integer",NT_l="integer",n_l="integer",phi_alpha="integer"), body = V2_CPP_rho_alpha.body,plugin="RcppArmadillo")



V2_CPP_log.diwish.partial.body<-'
RNGScope scope;
Rcpp::List Matrix_c(Matrix);
Rcpp::List Data_clist(Data_list);
NumericVector log_det=Matrix_log_det_input;
NumericVector dd=dof;
IntegerVector n_l = n;
int nn=n_l.size();
IntegerVector k_ll = k_l;
int K=k_ll.size();

NumericMatrix rr(nn,K);
NumericVector rtemp(1);

for (int k=0; k<K; k++){
SEXP matrix =  Matrix_c[k];
arma :: mat matrix_cc=Rcpp :: as < arma :: mat >(matrix);

for(int i=0; i<nn; i++) {
SEXP data_inv =  Data_clist[i];
arma :: mat data_inv_c=Rcpp :: as < arma :: mat >(data_inv);
arma :: mat op = data_inv_c * matrix_cc ;
rtemp[0]=arma::trace(op);
rr(i,k)=-100+0.5*dd[0]*log_det[k]- rtemp[0]*0.5;
}
}


NumericMatrix rr_output(clone(rr));
return(wrap(rr_output));


'


V2_CPP_log.diwish.partial<-cxxfunction(signature(Data_list="List",Matrix_log_det_input="numeric",Matrix="List",dof="numeric",k_l="integer",n="integer"), body = V2_CPP_log.diwish.partial.body,plugin="RcppArmadillo")





CPP_update.subject.level.body<-'
RNGScope scope;

IntegerVector G_subject = ll_s;  
IntegerVector Kll = K_l;  
NumericVector g_sum_rho = gg_rho_input;  
NumericVector pos = pos_wish;
NumericVector eee = eta_k_input;
Rcpp::List alist(adj_list_input);
int n = G_subject.size();
int K = Kll.size();
double PHIB = Rcpp::as<double>(phi_beta);

NumericVector prior(K);
NumericVector U(1);
NumericVector CDF(K);


for(int i=0; i<n; i++) {



SEXP ll = alist[i];
Rcpp::IntegerVector nei(ll);
int m=nei.size();  
for(int k=0; k<K; k++){
prior[k]=-100;
CDF[k]=0;
}
for(int j=0; j<m; j++){
prior[G_subject[nei[j]-1]-1]=prior[G_subject[nei[j]-1]-1]+PHIB;
}



CDF[0]=exp(prior[0]+g_sum_rho[i-1]+pos[i-1]+eee[0]);



for(int k=1; k<K; k++){

CDF[k]=CDF[k-1]+exp(prior[k]+g_sum_rho[k*n-1+i]+pos[k*n-1+i]+eee[k]);

}



U = CDF[K-1]*runif(1,0,1);
G_subject[i]=1;
for(int k=1; k<K; k++){
if(U[0]>CDF[k-1]){G_subject[i]=k+1;}
}



}


IntegerVector G_subject_output(clone(G_subject));
return(wrap(G_subject_output));


'

CPP_update.subject.level<-cxxfunction(signature(ll_s="integer",pos_wish="numeric" ,gg_rho_input="numeric",adj_list_input="List",K_l="integer",phi_beta="numeric",eta_k_input="numeric"),body = CPP_update.subject.level.body,plugin="RcppArmadillo")




CPP_group.likelihood.body<-'
RNGScope scope;

IntegerVector G_group = ll_g;  
IntegerVector Kll = K_l;  
NumericVector s_sum_rho = ss_rho_input;  
Rcpp::List alist(adj_list_input);
int n = G_group.size();
int K = Kll.size();
double PHIB = Rcpp::as<double>(phi_beta);

NumericVector prior(K);
NumericVector U(1);
NumericVector CC(K);


for(int i=0; i<n; i++) {



SEXP ll = alist[i];
Rcpp::IntegerVector nei(ll);
int m=nei.size();  
for(int k=0; k<K; k++){
prior[k]=-100;
CC[k]=0;
}
for(int j=0; j<m; j++){
prior[G_group[nei[j]-1]-1]=prior[G_group[nei[j]-1]-1]+PHIB;
}



CC[0]=prior[0]+s_sum_rho[i-1];



for(int k=1; k<K; k++){

CC[k]=CC[k-1]+prior[k]+s_sum_rho[k*n-1+i];

}

U[0]=U[0]+CC[K-1];


}
U[0]=U[0];


NumericVector U_output(clone(U));
return(wrap(U_output));


'

CPP_group.likelihood<-cxxfunction(signature(ll_g="integer",ss_rho_input="numeric" ,adj_list_input="List",K_l="integer",phi_beta="numeric"), body = CPP_group.likelihood.body,plugin="RcppArmadillo")



CPP_subject.likelihood.body<-'
RNGScope scope;

IntegerVector G_subject = ll_s;  
IntegerVector Kll = K_l;  
NumericVector g_sum_rho = gg_rho_input;  
NumericVector pos = pos_wish;
Rcpp::List alist(adj_list_input);
int n = G_subject.size();
int K = Kll.size();
double PHIB = Rcpp::as<double>(phi_beta);

NumericVector prior(K);
NumericVector CC(K);
NumericVector U(1);


for(int i=0; i<n; i++) {



SEXP ll = alist[i];
Rcpp::IntegerVector nei(ll);
int m=nei.size();  
for(int k=0; k<K; k++){
prior[k]=-100;
CC[k]=0;
}
for(int j=0; j<m; j++){
prior[G_subject[nei[j]-1]-1]=prior[G_subject[nei[j]-1]-1]+PHIB;
}



CC[0]=prior[0]+g_sum_rho[i-1]+pos[i-1];



for(int k=1; k<K; k++){

CC[k]=CC[k-1]+prior[k]+g_sum_rho[k*n-1+i]+pos[k*n-1+i];

}

U[0]=U[0]+CC[K-1];
}

U[0]=U[0];


NumericVector U_output(clone(U));
return(wrap(U_output));


'


CPP_subject.likelihood<-cxxfunction(signature(ll_s="integer",pos_wish="numeric" ,gg_rho_input="numeric",adj_list_input="List",K_l="integer",phi_beta="numeric"),body = CPP_subject.likelihood.body,plugin="RcppArmadillo")




######The main_function#######
Potts_Bayesian_Semi_Multi<-function(Data,net,n,p=3,N,
                                    K,T_index,iter,
                                    LU_alpha=c(0.25,3), 
                                    LU_beta=c(0.25,3), 
                                    LU_xi=c(0,1),
                                    rho_alpha=0.3,
                                    rho_beta=2,
                                    rho_xi=0.1){
  
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  adj_list<- apply(net==1,1,which)
  S=diag(3)
  
  eta_k=sapply(1:K,function(k) -k^rho_xi)
  
  ###label_sub and label_group is a label index for each subject and group, respectively
  largest_eigen<-function(m){log(det(m))}
  stat.Data_matrix<-sapply(1:N, function(NN) sapply(Data[[NN]],largest_eigen) )
  result<-Mclust(rowMeans(stat.Data_matrix),verbose =FALSE)
  class<-result$classification
  label_sub<-matrix(rep(class,N),n,N)
  label_group<-matrix(rep(class,2),n,2)
  remove(stat.Data_matrix,result,class)
  
  
  
  Bookkeeping=list()
  Bookkeeping$label_sub[[1]]<-label_sub
  Bookkeeping$label_group[[1]]<-label_group
  ###Matrix_list is a list containing the scale matrix of each cluster (including ones which are not existed in current state)
  Matrix_list<-lapply(1:K, function(k) rwish(10, diag(3)))
  Bookkeeping$Matrix_list[[1]]<-Matrix_list
  ###dof_iw is the dof of inverse Wishart
  dof_iw=5
  Bookkeeping$dof_iw[[1]]<-dof_iw
  ###dof_w is the dof of inverse Wishart
  dof_w=3
  Bookkeeping$dof_w[[1]]<-dof_w
  
  
  ori_inv.Data<-lapply(1:N, function(NN) list.apply(Data[[NN]],solve) )
  
  
  
  A_prior=diag(3)
  
  it=1
  
  
  
  for (it in 1:iter){
    
    ####The codes is written in the classic Wishart/invWishart. Thus we do some algebra modifications.
    
    inv.Data=lapply(1:N, function(NN) lapply(1:n,function(nn) ori_inv.Data[[NN]][[nn]]*(dof_iw-p-1)) )
    logdet.Data<-lapply(1:N, function(NN) sapply(1:n, function(nn) logdet(Data[[NN]][[nn]]/(dof_iw - p - 1)) ) )
    A_prior=diag(3)/dof_w
    
    
    
    
    choose=sample(c(0,1),1,TRUE,c(0.5,0.5))
    if(choose==1){
      gc(FALSE)
    }
    #Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, it)
    #######################################
    ###### Update Matrix List (Gibbs)######
    #######################################
    ## new matrix 
    #unique_label=1:K
    #matrix_new=lapply(lapply(unique_label,function(kk) Reduce("+", lapply(1:N, function(NN) list.sum(NN,kk)  ))+A_prior), solve)
    ## new degree of freedom
    #dof_new=sapply(unique_label, function(kk) dof_w+(dof_iw)*sum(label_sub==kk))
    
    
    matrix_new<-lapply(1:K, function(kk) A_prior)
    dof_new<-rep(dof_w,K)
    ## new matrix
    unique_label=sort(unique(as.vector(label_sub)))
    #unique_label=1:K
    matrix_new_temp<-lapply(1:length(unique_label), function(kk) A_prior)
    for (NN in 1:N){
      lls_input<-label_sub[,NN]
      inv.Data_input<-inv.Data[[NN]]
      matrix_new_temp<-CPP_new.matrix(lls_input,unique_label,inv.Data_input,matrix_new_temp)
      remove(lls_input,inv.Data_input)
    }
    tobeupdated<-lapply(matrix_new_temp,solve)
    matrix_new[unique_label]<-tobeupdated
    
    ## new degree of freedom
    
    lls_input<-(label_sub)
    dof_new.temp<-CPP_new.dof(lls_input,unique_label)
    dof_new[unique_label]<-dof_new.temp*dof_iw+dof_w
    remove(lls_input)
    
    ##Sample new matrix list
    Matrix_list<-lapply(1:K, function(kk) rwish(dof_new[[kk]],matrix_new[[kk]]))
    Matrix_list_log_det<-sapply(1:K, function(kk) log(det(Matrix_list[[kk]])) )
    
    
    ####################################
    ###### Update dof_w (MH)###########
    ####################################
    dof_w_can=rnorm(1,dof_w,1)
    if (dof_w_can<(p)){dof_w_can=p}
    if (dof_w_can>(20)){dof_w_can=20}

    can= sum(sapply(1:K, function(kk) dWishart(Matrix_list[[kk]], A_prior, dof_w_can,TRUE)))
    old= sum(sapply(1:K, function(kk) dWishart(Matrix_list[[kk]], A_prior, dof_w    ,TRUE)))

    prob=min(1,exp(can-old))

    if(is.infinite(prob)){
      if(can>old){
        prob=1
      }else{
        prob=0
      }
    }
    Accept=sample(c(0,1),1,prob = c(1-prob,prob))
    if(Accept==1){dof_w=dof_w_can}


    ####################################
    ###### Update dof_iw (MH)###########
    ####################################
    dof_iw_can=rnorm(1,dof_iw,1)
    if (dof_iw_can<(p+1)){dof_iw_can=p+2}
    if (dof_iw_can>20){dof_iw_can=20}

    #can= sum(sapply(1:N, function(NN) sapply(1:n, function(nn) dInverseWishart(Data[[NN]][[nn]],solve(Matrix_list[[label_sub[nn,NN]]]),dof_iw_can,TRUE) )) )
    #old= sum(sapply(1:N, function(NN) sapply(1:n, function(nn) dInverseWishart(Data[[NN]][[nn]],solve(Matrix_list[[label_sub[nn,NN]]]),dof_iw    ,TRUE) )) )


    can<-old<-0
    for(NN in 1:N){
      lmg_can<-lmgamma(0.5*dof_iw_can,p)
      added_can<-CPP_dof_iw(dof_iw_can,p,Matrix_list_log_det,logdet.Data[[NN]],label_sub[,NN],lmg_can,1:n)
      can<-can+added_can

      lmg<-lmgamma(0.5*dof_iw,p)
      added_old<-CPP_dof_iw(dof_iw    ,p,Matrix_list_log_det,logdet.Data[[NN]],label_sub[,NN],lmg,1:n)
      old<-old+added_old
    }

    prob=min(1,exp(can-old))

    if(is.infinite(prob)){
      if(can>old){
        prob=1
      }else{
        prob=0
      }
    }
    try(Accept<-sample(c(0,1),1,TRUE,c(1-prob,prob)),TRUE)


    if(Accept==1){dof_iw=dof_iw_can}
    
    
    
    ####################################
    ###### Update indicator label#######
    ####################################
    
    # ss_rho_all=matrix(0,n,K)
    # for(k in 1:K){
    #   #ss_rho_all[,k]<-sapply(1:n, function(nn)  sum(label_sub[nn,]==k))*rho_alpha
    #   llb<-(t(label_sub))
    #   try(ss_rho_all[,k]<-CPP_rho_alpha(llb,k,1:N,1:n,1)*rho_alpha,TRUE)
    # }
    
    llb<-(t(label_sub))
    ss_rho_all<-V2_CPP_rho_alpha(llb,1:K,1:N,1:n,rho_alpha)
    
    ss_rho_all_input<-(ss_rho_all)
    LG1_in<-label_group[,1]
    LG2_in<-label_group[,2]
    
    Diff<-CPP_indicator(ss_rho_all_input, LG1_in,  LG2_in,adj_list,1:n,1:K,rho_beta)
    #try(Ratio<-CPP_ratio(ss_rho_all,label_group[,1], label_group[,2],adj_list,1:n,1:K,rho_beta),TRUE)
    remove(llb, ss_rho_all,ss_rho_all_input,LG1_in,LG2_in)
    
    
    ####################################
    ###### Update Group Label###########
    ####################################
    # for (t in 1:2){
    #   for(nn in 1:n){
    #     label_group[n,t]<-sample(1:K,1,prob=exp(sapply(1:K, function(kk) update.group.level(nn,kk))))
    #   }
    # }
    
    for (t in 1:2){
      # ss_rho=matrix(0,n,K)
      # for(k in 1:K){   
      #   #ss_rho[,k]<-sapply(1:n, function(nn)  sum(label_sub[nn,which(T_index==t)]==k))*rho_alpha
      #   llb<-(t(label_sub[,which(T_index==t)]))
      #   try(ss_rho[,k]<-CPP_rho_alpha(llb,k,1:(N/2),1:n,1)*rho_alpha,TRUE)
      # }
      
      llb_input<-(t(label_sub[,which(T_index==t)]))
      ss_rho<-V2_CPP_rho_alpha(llb_input,1:K,1:(N/2),1:n,rho_alpha)
      ss_rho_input<-(ss_rho)
      LG_input<-label_group[,t]
      LG_output<-CPP_update.group.level(LG_input,ss_rho_input,adj_list,1:K,rho_beta)
      label_group[,t]<-LG_output
      
      remove(llb_input, ss_rho,ss_rho_input,LG_input,LG_output)
    }
    
    
    
    ####################################
    ###### Update Subject Label#########
    ####################################
    #pos_wish=rep(0,K)
    #pos_wish[1:K]=sapply(1:K, function(kk) dWishart(Matrix_list[[unique_label[kk]]],matrix_new[[kk]],dof_new[[kk]],TRUE))-100
    #pos_wish[1:K]<-0
    #pos_wish[1:K]=sapply(1:K, function(kk) sum(sapply(1:N, function(NN) log.diwish.partial(inv.Data[[NN]][[nn]],dof_iw,Matrix_list[[kk]]))))
    #pos_wish[setdiff(1:K,unique_label)]=sapply(setdiff(1:K,unique_label), function(kk) dWishart(Matrix_list[[kk]],A_prior,dof_w,TRUE))
    
    
    for (sub in 1:N){
      
      #pos_wish<-(sapply(1:K, function(kk) sapply(1:n, function(nn) dWishart(Data[[sub]][[nn]],Matrix_list[[kk]],dof_iw,TRUE))))
      
      #pos_wish<-sapply(1:K, function(kk) CPP_log.diwish.partial(inv.Data[[sub]],Matrix_list_log_det[kk],Matrix_list[kk],dof_iw,1:n))
      #pos_wish<-(pos_wish)
      
      #pos_wish<-matrix(0,n,K)
      #gg_rho<-matrix(0,n,K)
      # for (kk in 1:K){
      #   try(pos_wish[,kk]<-CPP_log.diwish.partial(inv.Data[[sub]],Matrix_list_log_det[kk],Matrix_list[kk],dof_iw,1:n),TRUE)
      #   #gg_rho[,kk]<-sapply(1:n, function(nn)  sum(label_group[nn,T_index[sub]]==kk))*rho_alpha
      #   try(gg_rho[,kk]<-CPP_rho_alpha(label_group[,T_index[sub]],k,1,1:n,1)*rho_alpha,TRUE)
      # }
      inv_Data_sub_input<-inv.Data[[sub]]
      pos_wish<-V2_CPP_log.diwish.partial(inv_Data_sub_input,Matrix_list_log_det,Matrix_list,dof_iw,1:K,1:n)
      gg_rho<-V2_CPP_rho_alpha(label_group[,T_index[sub]],1:K,1,1:n,rho_alpha)
      
      gg_rho_in<-(gg_rho)
      pos_wish_input<-(pos_wish)
      label_sub_in<-label_sub[,sub]
      label_sub[,sub]<-CPP_update.subject.level(label_sub_in,pos_wish_input,gg_rho_in,adj_list,1:K,rho_beta,eta_k)
      
      remove(pos_wish,gg_rho,pos_wish_input,gg_rho_in,inv_Data_sub_input,label_sub_in)
      
    }
    
    
    
    
    ####################################
    ######Exchange Algorithm############
    ####################################
    
    ####Step 1: generate candidate value
    
    
    for (jj in 1:3){
    
    u <- rnorm(3)*10^(-10)
    theta_prop <- log(c(rho_alpha,rho_beta,rho_xi)) + S %*% u
    
    
    
    
    rho_alpha_can<-exp(theta_prop[1])
    rho_beta_can<-exp(theta_prop[2])
    rho_xi_can<-exp(theta_prop[3])
    
    
    if(
      LU_alpha[1]<rho_alpha_can & rho_alpha_can<LU_alpha[2] &
       LU_beta[1]<rho_beta_can & rho_beta_can<LU_beta[2] &
       LU_xi[1]<rho_xi_can & rho_xi_can<LU_xi[2]){
      #if( lower<rho_alpha_can & rho_alpha_can<Ca & lower<rho_beta_can & rho_beta_can<Cb){
      
      eta_k_can=sapply(1:K,function(k) -k^rho_xi_can)
      
      
      
      ##Group##
      Y_label_group=label_group
      for (t in 1:2){
        # ss_rho=matrix(0,n,K)
        # for(k in 1:K){
        #   #ss_rho[,k]<-sapply(1:n, function(nn)  sum(label_sub[nn,which(T_index==t)]==k))*rho_alpha
        #   llb<-(t(label_sub[,which(T_index==t)]))
        #   try(ss_rho[,k]<-CPP_rho_alpha(llb,k,1:(N/2),1:n,1)*rho_alpha,TRUE)
        # }
        
        llb_input<-(t(label_sub[,which(T_index==t)]))
        ss_rho<-V2_CPP_rho_alpha(llb_input,1:K,1:(N/2),1:n,rho_alpha_can)
        
        ss_rho_input<-(ss_rho)
        Y_label_group[,t]<-CPP_update.group.level(label_group[,t],ss_rho_input,adj_list,1:K,rho_beta_can)
        
        remove(llb_input,ss_rho,ss_rho_input)
        
      }
      
      
      ##subject##
      Y_label_sub=label_sub
      for (sub in 1:N){
        
        #pos_wish<-(sapply(1:K, function(kk) sapply(1:n, function(nn) dWishart(Data[[sub]][[nn]],Matrix_list[[kk]],dof_iw,TRUE))))
        
        #pos_wish<-sapply(1:K, function(kk) CPP_log.diwish.partial(inv.Data[[sub]],Matrix_list_log_det[kk],Matrix_list[kk],dof_iw,1:n))
        #pos_wish<-(pos_wish)
        
        #pos_wish<-matrix(0,n,K)
        #gg_rho<-matrix(0,n,K)
        # for (kk in 1:K){
        #   try(pos_wish[,kk]<-CPP_log.diwish.partial(inv.Data[[sub]],Matrix_list_log_det[kk],Matrix_list[kk],dof_iw,1:n),TRUE)
        #   #gg_rho[,kk]<-sapply(1:n, function(nn)  sum(label_group[nn,T_index[sub]]==kk))*rho_alpha
        #   try(gg_rho[,kk]<-CPP_rho_alpha(label_group[,T_index[sub]],k,1,1:n,1)*rho_alpha,TRUE)
        # }
        inv_Data_sub_input<-inv.Data[[sub]]
        pos_wish<-V2_CPP_log.diwish.partial(inv_Data_sub_input,Matrix_list_log_det,Matrix_list,dof_iw,1:K,1:n)
        gg_rho<-V2_CPP_rho_alpha(label_group[,T_index[sub]],1:K,1,1:n,rho_alpha_can)
        
        gg_rho_in<-(gg_rho)
        pos_wish_input<-(pos_wish)
        label_sub_in<-label_sub[,sub]
        Y_label_sub[,sub]<-CPP_update.subject.level(label_sub_in,pos_wish_input,gg_rho_in,adj_list,1:K,rho_beta_can,eta_k_can)
        
        remove(pos_wish,gg_rho,gg_rho_in,pos_wish_input,inv_Data_sub_input,label_sub_in)
      }
      
      ######Step 3: Acceptance ratio
      ####Compute group likelihood
      ####X
      ##old
      X_group_old<-X_group_can<-0
      for(t in 1:2){
        llb_input<-(t(label_sub[,which(T_index==t)]))
        ss_rho_all<-V2_CPP_rho_alpha(llb_input,1:K,1:(N/2),1:n,rho_alpha)
        ss_rho_all_input<-(ss_rho_all)
        added_old<-CPP_group.likelihood(label_group[,t],ss_rho_all_input,adj_list,1:K,rho_beta)
        X_group_old<-X_group_old+added_old
        remove(llb_input,ss_rho_all,ss_rho_all_input)
        
        ##can
        llb_input<-(t(label_sub[,which(T_index==t)]))
        ss_rho_all<-V2_CPP_rho_alpha(llb_input,1:K,1:(N/2),1:n,rho_alpha_can)
        ss_rho_all_input<-(ss_rho_all)
        added_can<-CPP_group.likelihood(label_group[,t],ss_rho_all_input,adj_list,1:K,rho_beta_can)
        added_can=added_old/rho_alpha*rho_alpha_can
        X_group_can<-X_group_can+added_can
        remove(llb_input,ss_rho_all,ss_rho_all_input)
      }
      ####Y
      ##old
      Y_group_old<-Y_group_can<-0
      for(t in 1:2){
        llb_input<-(t(Y_label_sub[,which(T_index==t)]))
        ss_rho_all<-V2_CPP_rho_alpha(llb_input,1:K,1:(N/2),1:n,rho_alpha)
        ss_rho_all_input<-(ss_rho_all)
        added_old<-CPP_group.likelihood(Y_label_group[,t],ss_rho_all_input,adj_list,1:K,rho_beta)
        Y_group_old<-Y_group_old+added_old
        remove(llb_input,ss_rho_all,ss_rho_all_input)
        
        ##can
        llb_input<-(t(Y_label_sub[,which(T_index==t)]))
        ss_rho_all<-V2_CPP_rho_alpha(llb_input,1:K,1:(N/2),1:n,rho_alpha_can)
        ss_rho_all_input<-(ss_rho_all)
        added_can<-CPP_group.likelihood(Y_label_group[,t],ss_rho_all_input,adj_list,1:K,rho_beta_can)
        Y_group_can<-Y_group_can+added_can
        remove(llb_input,ss_rho_all,ss_rho_all_input)
      }
      
      
      #####Compute subject likelihood
      ####X
      ##old
      X_sub_old<-0
      pos_wish_input<-(matrix(0,n,K))
      for (sub in 1:N){
        gg_rho<-V2_CPP_rho_alpha(label_group[,T_index[sub]],1:K,1,1:n,rho_alpha)
        gg_rho_in<-(gg_rho)
        added<-CPP_subject.likelihood(label_sub[,sub],pos_wish_input,gg_rho_in,adj_list,1:K,rho_beta)
        X_sub_old<-X_sub_old+added
        
        
      }
      
      remove(gg_rho_in,pos_wish_input)
      ##can
      X_sub_can<-0
      pos_wish_input<-(matrix(0,n,K))
      for (sub in 1:N){
        gg_rho<-V2_CPP_rho_alpha(label_group[,T_index[sub]],1:K,1,1:n,rho_alpha_can)
        gg_rho_in<-(gg_rho)
        added<-CPP_subject.likelihood(label_sub[,sub],pos_wish_input,gg_rho_in,adj_list,1:K,rho_beta_can)
        X_sub_can<-X_sub_can+added
        
      }
      remove(gg_rho_in,pos_wish_input)
      ####Y
      ##old
      Y_sub_old<-0
      pos_wish_input<-(matrix(0,n,K))
      for (sub in 1:N){
        gg_rho<-V2_CPP_rho_alpha(Y_label_group[,T_index[sub]],1:K,1,1:n,rho_alpha)
        gg_rho_in<-(gg_rho)
        added<-CPP_subject.likelihood(Y_label_sub[,sub],pos_wish_input,gg_rho_in,adj_list,1:K,rho_beta)
        Y_sub_old<-Y_sub_old+added
        
        
      }
      remove(gg_rho_in,pos_wish_input)
      ##can
      Y_sub_can<-0
      pos_wish_input<-(matrix(0,n,K))
      for (sub in 1:N){
        gg_rho<-V2_CPP_rho_alpha(Y_label_group[,T_index[sub]],1:K,1,1:n,rho_alpha_can)
        gg_rho_in<-(gg_rho)
        added<-CPP_subject.likelihood(Y_label_sub[,sub],pos_wish_input,gg_rho_in,adj_list,1:K,rho_beta_can)
        Y_sub_can<-Y_sub_can+added
        
      }
      remove(gg_rho_in,pos_wish_input)
      
      ratio=(X_sub_can+X_group_can+Y_sub_old+Y_group_old)-(Y_sub_can+Y_group_can+X_sub_old+X_group_old)
      
      prob=min(1,exp(ratio))
    
      if(is.infinite(prob)){
        if(ratio>0){
          prob=1
        }else{
          prob=0
        }
      }
      
      Accept<-sample(c(0,1),1,TRUE,c(1-prob,prob))
      
      
      
      
      ####################################
      ####################################
      ####################################
      
      # ss_rho_all=matrix(0,n,K)
      # for(k in 1:K){
      #   #ss_rho_all[,k]<-sapply(1:n, function(nn)  sum(label_sub[nn,]==k))*rho_alpha
      #   llb<-(t(label_sub))
      #   try(ss_rho_all[,k]<-CPP_rho_alpha(llb,k,1:N,1:n,1)*rho_alpha,TRUE)
      # }
      
      llb<-(t(label_sub))
      ss_rho_all<-V2_CPP_rho_alpha(llb,1:K,1:N,1:n,rho_alpha)
      
      ss_rho_all_input<-(ss_rho_all)
      LG1_in<-label_group[,1]
      LG2_in<-label_group[,2]
      
      Diff_temp<-CPP_indicator(ss_rho_all_input, LG1_in,  LG2_in,adj_list,1:n,1:K,rho_beta)
      #try(Ratio<-CPP_ratio(ss_rho_all,label_group[,1], label_group[,2],adj_list,1:n,1:K,rho_beta),TRUE)
      remove(llb, ss_rho_all,ss_rho_all_input,LG1_in,LG2_in)
      ####################################
      ####################################
      ####################################
      
      
      
      
      
      
      if(mean(Diff_temp)>0.4){
        Accept=0
        prob=0
      }
      
      #if(it<2000){
      #S<- ramcmc::adapt_S(S, u, prob, it )
      #}
      
      if(Accept==1 ){
        rho_alpha=rho_alpha_can
        rho_beta=rho_beta_can
        rho_xi=rho_xi_can
        eta_k=eta_k_can
      }
      
      
    }
    
    
    }
    
    
    
    if(it==1){
      Diff_all<-list.rbind(lapply(1:iter,function(it) Diff))
      rho_all<-list.rbind(lapply(1:iter,function(it) c(rho_alpha,rho_beta,rho_xi)))
    }
    
    
    #Bookkeeping$label_sub[[it]]<-label_sub
    #Bookkeeping$label_group[[it]]<-label_group
    #Bookkeeping$Diff[[it]]<-Diff
    #Bookkeeping$rho[[it]]<-c(rho_alpha,rho_beta,rho_xi)
    #Bookkeeping$Matrix_list[[it]]<-Matrix_list    
    #Bookkeeping$dof_iw[[it]]<-dof_iw
    #Bookkeeping$dof_w[[it]]<-dof_w
    Diff_all[it,]<-Diff
    rho_all[it,]<-c(rho_alpha,rho_beta,rho_xi)
    
    
    
     if(it%%10==0){
    #   print(it)
    #   print(sum(Diff))
    #   print(dof_iw)
    #   print(dof_w)
   #    print(c(rho_alpha,rho_beta,rho_xi))
    image(matrix(Diff,c,c))
    #   # choose=sample(c(0,1),1,TRUE,c(0.5,0.5))
    #   # if(choose==1){
    #   #   
    #   # }
    #   #
    #   #source('/Volumes/Transcend/Bayesian_Nonparametric_DTI_Simulations/functions_Hierachical_Potts_Bayes.R')
    #   #save.image("/mnt/home/zlan/Bayes_DTI/Ratio_real_data_K100.RData")
   }
    
    
  }
  close(pb)
  
  
  Bookkeeping$Diff<-Diff_all
  Bookkeeping$rho<-rho_all
  result_in<-duplicate(Bookkeeping, shallow = FALSE)
  return(result_in)
  
}











