
# Defining functions 



# Joint maximum likelihood estimation 

ipf_JML <- function(data,parinit=NA, 
                    iter=100, 
                    epsilon=0.000001) {
  # Data prep 
  
  n_pers <- nrow(data)

  # observed table margins 
  V <- rowSums(data) # row sums 
  W <- colSums(data) # columns sums

  n_item <- ncol(data)# number of columns 
  
  #Computing the expected sufficient marginals. 
  p_ij <- function(xi,delta) { 
    xi*delta/(1+xi*delta)
  }
  r_pij <-function(xi,delta,f){
    n <- length(xi)
    vec <- numeric(n)
    for (i in 1:n){
      vec[i] <- sum(f(xi[i],delta))
    }
    vec
  }
  
  c_pij <- function(xi,delta,f){
    n <- length(delta)
    vec <- numeric(n)
    for (i in 1:n){
      vec[i] <- sum(f(xi,delta[i]))
    }
    vec
  }
  
  
  # Item itfficulty parameters sum to 1. 
  prod1 <- function(par) {
  if(prod(par)==1){cbind(par,1)}
  else {
    geometric_mean <- exp(mean(log(par)))
    par <- par / geometric_mean
    cbind(par,geometric_mean)
    
  }
}
  
  # Defining initial parameters. 
  if (is.na(parinit)){
    xi <- rep(1,n_pers)
    delta <- rep(1,n_item)
    }
  
  else {
    xi <- parinit[,1]
    delta <- parinit[,2]
  }

  
  R <- r_pij(xi,delta,p_ij)
  C <- c_pij(xi,delta,p_ij)
  temp <- prod1(delta)

  
  mat_fun <- function(xi,delta,tem){
    mat <- matrix(ncol=n_item, nrow=n_pers)
  
    for ( i in 1:n_item){
      for ( j in 1:n_pers) {
        mat[j,i] <- p_ij(xi[j]*tem,delta[i])
      }
      
    }
    mat
    }
  
  # Resrict item difficulties to have product 1. 
 
  prod1 <- function(par) {
    if(prod(par)==1){cbind(par,1)}
    else {
      geometric_mean <- exp(mean(log(par)))
      par <- par / geometric_mean
      cbind(par,geometric_mean)
      
    }
  }

  for (t in 1:iter){
    #Update person parameter
    
    xi <- (V/R)*xi
   
     # Expected column sums 
    C <- c_pij(xi*temp[1,2],delta,p_ij)
    
    temp <- prod1(c(W/C)*delta)


    delta <- temp[,1]

    # Expected rowsums  
    R <- r_pij(xi*temp[1,2],delta,p_ij)
    

  } 
  
  list("Xi" = xi*temp[1,2],"Delta" = delta,
       "R" = R,"C" = C,
       "Iter" = mat_list)
  
}


# Joing maximum likelihood estimation where persons are grouped. 

JML_IFP <- function(data,parinit=NA, 
                    maxiter=1000,
                    epsilon=0.0001) {
  
  n_r <- table(rowSums(data))
  r <- as.integer(names(n_r)) # Score in each scoregroup 
  n_r <- unname(n_r)
  
  # target table margins 
  V <- n_r*r # row sums 
  
  W <- colSums(data) # columns sums

  n_item <- ncol(data) # number of columns 
  n_score <- length(n_r) 
  
  
  ## Observed sufficient margins 
  r_hat_k <- function(alpha,beta,f){
    n <- length(alpha)
    vec <- numeric(n)
    for (i in 1:n){
      vec[i] <- sum(beta*f(alpha[i],beta,n_r[i]))
    }
    alpha*vec
  }
  
  c_hat_k <- function(alpha,beta,f){
    n <- length(beta)
    vec <- numeric(n)
    for (i in 1:n){
      vec[i] <- sum(alpha*f(alpha,beta[i],n_r))
    }
    beta*vec
  }
  
  e_rj <- function(xi,delta,nr) {
    nr/(1+xi*delta)
  }

  # restriction of item difficulties
  prod1 <- function(par) {
    if(prod(par)==1){cbind(par,1)}
    else {
      geometric_mean <- exp(mean(log(par)))
      par <- par / geometric_mean
      cbind(par,geometric_mean)
      
    }
  }
  mat_fun <- function(xi,delta,tem){
    mat <- matrix(ncol=n_item, nrow=n_score)
    
    for ( i in 1:n_item){
      for ( j in 1:n_score) {
        mat[j,i] <- xi[j]*tem*delta[i]*e_rj(xi[j]*tem,delta[i],nr=n_r[j])
      }
      
    }
    
    mat <- cbind(mat,rowSums(mat))
    mat <- rbind(mat,colSums(mat))
    mat
  }
  
  # Running algorithm 
  if (is.na(parinit)){
    xi <- rep(1,n_score)
    delta <- rep(1,n_item)
  }
  
  else {
    xi <- parinit[,1]
    delta <- parinit[,2]
  }
  
  R <- r_hat_k(xi,delta,e_rj)
  C <- c_hat_k(xi,delta,e_rj)
  temp <- prod1(delta)
  
  for (t in 1:maxiter){

    if(sum((R - V)^2) < epsilon & 
       sum((C - W)^2) < epsilon ){break}
    xi <- c(V/R)*xi
 
    C <- c_hat_k(xi*temp[1,2],delta, e_rj)
 
    
    temp <- prod1(c(W/C)*delta)
    delta <- temp[,1]
    
    R <- r_hat_k(xi*temp[1,2],delta,e_rj)
    
  } 
  list("delta" = -log(delta),
       "data" =data,
       "r"=r,
       "V"=V,
       "W"=W,
       "i"= t, 
       "eps" = epsilon, 
       "table" = mat_fun(xi,delta,temp[1,2]))
  
}

# Conditional maximum likelihood estimation 

CML_IPF <- function(data,parinit=NA, 
                    maxiter=1000, 
                    epsilon=0.0001) {
  

  n_r <- table(rowSums(data))
  r <- as.integer(names(n_r)) # Score in each scoregroup 
  n_r <- unname(n_r)
  
  # target table margins 
  
  W <- colSums(data) # columns sums
  V <- n_r*r 
  n_item <- ncol(data) # number of columns 
  n_score <- length(n_r) 

  # Restriction of item difficulties. 
  prod1 <- function(par) {
    if(prod(par)==1){cbind(par,1)}
    else {
      geometric_mean <- exp(mean(log(par)))
      par <- par / geometric_mean
      cbind(par,geometric_mean)
      
    }
  }

  # Gamma function
  sum_prod <- function(r,delta){
    temp <- combn(delta,r,FUN=prod)
    sum(temp)
  }
  
  gamma_fun <- function(r_p,delta,k){
    if(r_p==0){x<- 1}
    else if (r_p==1){x<- sum(delta)}
    else if(r_p==k){x<-1}
    else {x <- sum_prod(r_p,delta)}
    x
  }
  
  
  # expected column sum and entries
  e_ri_1 <- function(delta,nr,r,i,k=n_item) {
    nr*(gamma_fun(r-1,delta[-i],k)/gamma_fun(r,delta,k))
  }
  
  
  
  c_hat_k_1 <- function(beta,f){
    n <- length(n_r)
    n_beta <- length(beta)
    vec <- numeric(n_beta)
    for (j in 1:n_beta){
      for (i in 1:n){
        vec[j] <- vec[j]+e_ri_1(beta,n_r[i],r[i],j)
        
      }
    }
    beta*vec
  }
  
  
  # Running algorithm 
  

  delta <- rep(1,n_item) # initial parameters 
  C <- c_hat_k_1(delta,e_ri_1)

  # expected table 
  mat_fun <- function(delta){
    mat <- matrix(ncol=n_item, nrow=n_score)
    
    for ( i in 1:n_item){
      for ( j in 1:n_score) {
        mat[j,i] <- delta[i]*e_ri_1(delta,n_r[j],r[j],i)
      }
      
    }
    mat <- cbind(mat,rowSums(mat))
    mat <- rbind(mat,colSums(mat))
    mat
  }
  
  for (t in 1:maxiter){
    if (is.na(sum((C - W)^2) < epsilon ))
      {warning("Model failure")
       delta <- rep(NA,n_item)
       break    }
    if(sum((C - W)^2) < epsilon ){
  
      break
      }
    
    delta <- prod1(W/C*delta)[,1]
    C <- c_hat_k_1(delta,e_ri_1)
  }
  list("delta" = -log(delta),
       "r" = r,
       "V" = V,
       "W" = W,
       "data" =data,
       "i"= t, 
       "eps" = epsilon, 
       "table" = mat_fun(delta))
  
}



# Wrapper function 
Rasch_IPF<- function(data,
                     method="CML", 
                     parinit=NA, 
                     maxiter=1000, 
                     removeExtreme = TRUE,
                     epsilon=0.0001){
  
  if(removeExtreme){data <- data[!(rowSums(data)==0|rowSums(data)==max(rowSums(data))),]}
  if(any(rowSums(data)==0)){stop("If removeExtreme = FALSE extreme data points must be removed")}
  
  if(method=="CML"){ 
  result <- CML_IPF(data,
                    parinit, 
                    maxiter, 
                    epsilon)
  }
  if(method =="JML"){
    result<- JML_IFP(data,
                     parinit, 
                     maxiter, 
                     epsilon)
  }
  if(!(method=="JML"|method=="CML")){stop("Undefined method")}
  result$method <- method
  
  structure(result,class="Rasch_IPF_Estimation")
  
  
}

# Print funciton for dichotomous rasch model estimation 
print.Rasch_IPF_Estimation <- function(list){
  mat <- list$table
  dimnames(mat) <- list(c(list$r,"C"),c(names(list$delta),"R"))
  est <- t(t(c(list$delta)))
  colnames(est) <- " "
  cat(" \n",
      "Iterative proportional fitting Rasch model:\n Estimates and table \n
      Item difficulty parameters\n" )
  print(est)
  cat("\nTable:\n")
  print(mat)
  cat("\n Row target: ",list$V,
      "\n Column target: ", list$W)
  cat("\n \n Iterations:" ,list$i,"\n Method:", list$method)
}



## Certainty of parameters. 
Std_error_BS <- function(obj,data, B=1000){
  n <- nrow(data)
  n_item <-ncol(data)
  delta_mat <- matrix(ncol=n_item,nrow=B)
  for (i in 1:B){
    print(i)
    temp <- data[sample(1:n,n,replace=TRUE),]
    est <- Rasch_IPF(temp,method=obj$method, removeExtreme = TRUE, epsilon = obj$eps)
    delta_mat[i,] <- est$delta
  }
  delta_mat <- delta_mat
  sd <- apply(delta_mat, 2, function(x){sd(x,na.rm=TRUE)})
  structure(list("delta_sim" = delta_mat,
                 "delta" =obj$delta,
                 "sd" = sd,
                 "B" = B), class = "IPF_rasch_CI")
  }
                 
 



print.IPF_rasch_CI <- function(list){
  cat(" \n",
      "Iterative proportional fitting Rasch model:\n Estimates and confidence intervals with method bootstrap \n") 
      print(matrix(c(list$delta,
               list$sd,
               list$delta-1.96*list$sd,
               list$delta+1.96*list$sd),ncol=4,
             dimnames=list(names(list$delta),
                          c("Estimate", "Std. Error", "lower CI", "upper CI"))))
      
   cat("\n Repititions:" ,list$B)
}



autoplot.IPF_rasch_CI <- function(obj, grid=FALSE) {
  n <- length(obj$delta)
  d_plots <- list()
    for (i in 1:n){
      d_plots[[i]] <- ggplot() + geom_histogram(aes(x=obj$delta_sim[,i],y=..density..), col="black", fill="orange") + 
        geom_vline(aes(xintercept= obj$delta[i]),color="red") + 
        geom_vline(aes(xintercept= mean(obj$delta_sim[,i],na.rm=TRUE))) +
        geom_vline(aes(xintercept=quantile(obj$delta_sim[,i],0.025,na.rm=TRUE)), linetype="dashed")+ 
        geom_vline(aes(xintercept=quantile(obj$delta_sim[,i],0.975,na.rm=TRUE)), linetype="dashed") + theme_bw()
    
    }
  
  if(grid){do.call("grid.arrange", c(d_plots, ncol=2))}
  
  
  else{
    for ( i in 1:n){
      print(d_plots[[i]])
    }}
}
  

