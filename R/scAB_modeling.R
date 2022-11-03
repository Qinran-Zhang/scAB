###  Non-negative Matrix Factorization
#' Classical non-negative matrix factorization algorithm
#'
#' @param X  a non-negative matrix
#' @param K  the rank of matrix factorization
#' @param maxiter the maximum number of iterations
#'
#' @return a list with the submatrix and loss value
#' @export
#'
#' @examples
NMF <- function(X, K,maxiter=2000){
  eps=2.2204e-256
  nr = dim(X)[1];
  nc = dim(X)[2];
  W = matrix(runif(nr*K),nrow=nr, ncol=K);
  H = matrix(runif(K*nc),nrow=K, ncol=nc);
  loss_func=function(X,W,H){
    loss<-norm(X-W%*%H,"F")^2
    return(as.numeric(loss))
  }
  
  for (iter in 1:maxiter){
    H = H*(t(W)%*%X)/ ((t(W)%*%W)%*%H )
    W=W*(X%*%t(H))/(W %*% H %*% t(H)  )
    {
      if(iter!=1){
        eucl_dist = loss_func(X,W,H)
        d_eucl=abs(eucl_dist-old_eucl)
        if(d_eucl<10^(-5)) {break;}
        old_eucl=eucl_dist
      }
      else{ old_eucl =loss_func(X,W,H)
      }
    }
    iter=iter+1
  }
  return(list(W=W,H=H,iter=iter,loss=loss_func(X,W,H) ))
}

###  Non-negative Matrix Factorization with phenotype and cell-cell similarity regularization.
#' Non-negative Matrix Factorization with phenotype and cell-cell similarity regularization, for identifing phenotype-associated cell states at different resolutions.
#'
#' @param Object  a scAB_data object
#' @param K  the rank of matrix factorization
#' @param maxiter the maximum number of iterations
#' @param alpha Coefficient of phenotype regularization
#' @param alpha_2 Coefficient of cell-cell similarity regularization
#'
#' @return a list with the submatrix and loss value
#' @export
#'
#' @examples
scAB <- function(Object, K,alpha=0.005,alpha_2=0.005,maxiter=2000){
  seed = ifelse(Object$method=="survival",7,5)
  if(Object$method!="") set.seed(seed)
  X <- Object$X
  A <- Object$A
  L <- Object$L
  D <- Object$D
  S <- Object$S
  eps=2.2204e-256
  nr = dim(X)[1];
  nc = dim(X)[2];
  W = matrix(runif(nr*K),nrow=nr, ncol=K);
  H = matrix(runif(K*nc),nrow=K, ncol=nc);
  loss_func=function(X,W,H,S,L,alpha,alpha_2){
    loss<-norm(X-W%*%H,"F")^2 +
      alpha*( norm(S%*%W,"F")^2) +
      alpha_2*sum(diag(H %*% L %*% t(H)))
    return(as.numeric(loss))
  }
  tD <- t(D)
  tA <- t(A)
  for (iter in 1:maxiter){
    W_old <- W
    H <- H*(t(W)%*%X+alpha_2*H%*%tA)/ ((t(W)%*%W)%*%H + alpha_2*H%*%tD)
    Pena <- (S%*%S)%*%W
    W <- W*(X%*%t(H))/(W %*% H %*% t(H)  + alpha*Pena )
    {
      if(iter!=1){
        eucl_dist <- loss_func(X,W,H,S,L,alpha,alpha_2)
        d_eucl <- abs(eucl_dist-old_eucl)
        if(d_eucl<10^(-5)) {break;}
        old_eucl <- eucl_dist
      }
      else{ old_eucl <- loss_func(X,W,H,S,L,alpha,alpha_2)
      }
    }
    iter <- iter+1
  }
  return(list(W=W,H=H,iter=iter,loss=loss_func(X,W,H,S,L,alpha,alpha_2),method=Object$method ))
}


###  Create subsets of cross-validation
#'
#' @param k  k-fold cross validation
#' @param datasize  the size of samples
#' @param seed random seed
#'
#' @return a list with subsets of cross-validation
#' @export
#'
#' @examples
CVgroup <- function(k,datasize,seed = 0){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]
  temp <- sample(n,datasize)
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])
  return(cvlist)
}


###  Selection of parameter K
#'
#' @param Object a scAB_data object
#' @param k_max  the maximum value of the rank in the matrix factorization
#' @param repeat_times  the number of repetitions
#' @param seed random seed
#' @param verbose Print output
#'
#' @return an appropriate value of K
#' @export
#'
#' @examples
select_K <- function(Object, K_max=20, repeat_times=10, maxiter=2000, seed=0, verbose = FALSE){
  
  X=Object$X
  set.seed(seed)
  K_all=2:K_max
  dist_K=matrix(NA,nrow=length(K_all),ncol=repeat_times)
  dist_all=norm(X,"F")
  eii=c()
  for(Ki in K_all)
  {
    for(Kj in 1:repeat_times)
    {
      res_ij=NMF(X=X, K=Ki,maxiter=maxiter)
      dist_K[Ki-1,Kj]=norm(X-res_ij$W%*%res_ij$H,"F")^2
      if(Kj==repeat_times & verbose) print(paste0("loss of ",Ki,": ",mean(dist_K[Ki-1,])))
    }
    if(Ki==2) next;
    eii[Ki-1]= (rowMeans(dist_K)[Ki-2]-rowMeans(dist_K)[Ki-1])/(rowMeans(dist_K)[1]- rowMeans(dist_K)[Ki-1])
    if(rowMeans(dist_K)[Ki-2]-rowMeans(dist_K)[Ki-1]<=0) break;
    if(eii[Ki-1]<0.05 ) break;
  }
  K=Ki-1
  return(K)
}


###  Selection of parameter alpha and alpha_2
#'
#' @param Object a scAB_data object
#' @param k  the rank of matrix factorization
#' @param cross_k  k-fold cross validation
#' @param seed random seed
#'
#' @return a list with alpha, alpha_2 and cross-validation error
#' @importFrom MASS ginv
#' @import survival
#' @export
#'
#' @examples
select_alpha <- function(Object,K,cross_k=5,seed=0){
  if(Object$method=="survival")  train_phenotype <- Object$phenotype
    else {
      train_phenotype=data.frame(status=ifelse(Object$phenotype,1,0),time=ifelse(Object$phenotype,1,100))
      rownames(train_phenotype)<-rownames(Object$X)
    }
  train_data <- Object$X
  A_cv <- Object$A
  L_cv <- Object$L
  D_cv <- Object$D
  datasize <- nrow(train_data)
  cvlist <- CVgroup(k = cross_k,datasize = datasize,seed = seed)
  para_1_list <- c(0.01,0.005,0.001)
  para_2_list <- c(0.01,0.005,0.001)
  result_cv <- matrix(NA,nrow=length(para_1_list),ncol=length(para_2_list))
  times_para <- 0
  pb <- txtProgressBar(style=3)
  for(para_1 in 1:length(para_1_list)){
    for(para_2 in 1:length(para_2_list)){
      cv_c <- c()
      for(cvi in 1:cross_k)
      {
        train <- train_data[-cvlist[[cvi]],]
        test <- train_data[cvlist[[cvi]],]
        train_c <- train_phenotype[-cvlist[[cvi]],]
        test_c <- train_phenotype[cvlist[[cvi]],]
        train <- train/norm(train,"F")
        ss <- guanrank(train_c[,c("time","status")])
        S <- diag(1-ss[rownames(train_c),3])# *N
        Object_cv <- list(X=train,S=S,phenotype=train_c,A=A_cv,L=L_cv,D=D_cv,method="")
        class(Object_cv) <- "scAB_data"
        s_res <- scAB(Object=Object_cv, K=K,alpha=para_1_list[para_1],alpha_2=para_2_list[para_2],maxiter=2000)
        ginvH <- MASS::ginv(s_res$H)
        new_W <- test%*%ginvH
        clin_km <- data.frame(time=train_c$time,status=train_c$status,s_res$W)
        res.cox <- survival::coxph(survival::Surv(time, status) ~., data = clin_km)
        pre_test <- predict(res.cox,data.frame(new_W))
        cv_c[cvi] <- survival::concordance(survival::coxph(survival::Surv(test_c$time, test_c$status) ~pre_test))$concordance
      }
      result_cv[para_1,para_2] <- mean(cv_c)
      times_para <- times_para+1
      setTxtProgressBar(pb, times_para/length(para_1_list)/length(para_2_list))
    }
  }
  close(pb)
  
  para_index <- as.numeric(which(result_cv==max(result_cv),arr.ind=TRUE))
  alpha_1 <- para_1_list[para_index[1]]
  alpha_2 <- para_2_list[para_index[2]]
  return(para=list(alpha_1=alpha_1,alpha_2=alpha_2,result_cv=result_cv))
}


###  identification cells above the threshold
#'
#' @param W individual matrix
#' @param H cell matrix
#' @param tred threshold
#'
#' @return a list with cells above the threshold
#' @importFrom  diptest dip.test
#' @importFrom multimode locmodes
#'
#' @export
#'
#' @examples
findModule <- function(H,tred=2,do.dip=FALSE){
  K = dim(H)[1]
  I=length(H)
  module=list()
  meanH=rowMeans(H)
  sdH=apply(H,1,sd)
  for(i in 1:K){
    x <- H[i,]
    if( do.dip & diptest::dip.test(x)$p.value < 0.05)
    { modes <- multimode::locmodes(x, mod0 = 2)
      module=c(module,list(which( x > modes$locations[2])))}
    else {module=c(module,list( which( H[i,]-meanH[i] > tred*sdH[i] )) )}
  }
  return(module)
}



###  Subsets identification
#'
#' @param Object a Seurat object
#' @param scAB_Object a scAB_data object
#' @param tred threshold
#'
#' @return a Seurat object
#' @import Seurat
#' @export
#'
#' @examples
findSubset <- function(Object, scAB_Object, tred = 2){
  do.dip <- ifelse(scAB_Object$method=="binary",1,0)
  module <- findModule(scAB_Object$H, tred = tred, do.dip = do.dip)
  scAB_index <- unique(unlist(module))
  
  scAB_select <- rep("Other cells",ncol(Object))
  scAB_select[scAB_index] <- "scAB+ cells"
  Object <- Seurat::AddMetaData(Object,metadata = scAB_select, col.name = "scAB_select")
  
  for(i in 1:length(module)){
    M <- rep("Other cells",ncol(Object))
    M[as.numeric(module[[i]])]="scAB+ cells"
    Object <- Seurat::AddMetaData(Object,metadata = M, col.name = paste0("scAB_Subset",i))
    Object <- Seurat::AddMetaData(Object,metadata = scAB_Object$H[i,], col.name = paste0("Subset",i,"_loading"))
  }
  return(Object)
}

