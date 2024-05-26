CML_IPF_PCM <- function(data,
                        removeExtreme = TRUE,
                        x=NA,
                        parinit=NA, 
                    maxiter=1000, 
                    epsilon=0.0001 ,
                    return_tab = TRUE) {
  data <- as.matrix(data)
  
  if(is.na(x)){
    m <- max(data)
    x <- seq(1,m)} # if answer categories are not supplied they are found in data
  
  n_item <- ncol(data) # number of columns = items
  
  
  if(removeExtreme){data <- data[!(rowSums(data)==0|rowSums(data)==m*n_item),]} # remove observations where all replies are 0 or m 
  

  n_r <- table(rowSums(data))
  r <- as.integer(names(n_r)) # Score in each scoregroup 
  n_r <- unname(n_r)
  n_score <- length(n_r) 
  
  
  # Computing the observed sufficient marginals 
  W <- c()
  for (xhat in x){
    data_temp <- (data==xhat)*xhat
    W <- c(W,colSums(data_temp))
  }
  V <- n_r*r 

  #Gamma funcion is computed recursively. 
  gamma_rec <- memoise::memoise(function(beta,r){
    if(is.vector(beta)) {
      k <- length(beta)
      m <- 1
      beta <- as.matrix(beta)
    }
    else{k<- nrow(beta)
    m <- ncol(beta)}
    if (r==0){return(1)}
    if (r==1){return(sum(beta[,1]))}
    if (r==(k*m)){return(prod(beta[,m]))}
    if(r<0|nrow(beta)==0){return(0)}
    
    else{
      i <- 1
      temp <- 0 
      for (x in 1:m){
        temp <- temp + beta[i,x]*gamma_rec(matrix(beta[-i,],ncol=m),r-x)
      }
      return(temp + gamma_rec(matrix(beta[-i,],ncol=m),r))        
    }
    
  })
  
  # function to make sure the item difficulty parameters sum to 0. 
  prod1_pc <- function(matrix) {
    if(prod(matrix[,ncol(matrix)])==1){return(matrix)}
    else {
      geometric_mean <- exp(mean(log(matrix)))
      matrix <- (matrix / geometric_mean)
      
      
      return(matrix)
    }
    
  }
  
  
  
  
  update_beta<-function(X,beta){
    beta_upt <- X*as.vector(beta)
    prod1_pc(matrix(beta_upt,ncol=m))
  }
  
  
  #expected cell probabilities
  e_ri_1 <- function(beta,nr,r,i,x) {
     nr*x*beta[i,x]*(gamma_rec(beta[-i,],r-x)/gamma_rec(beta,r))
  }  
  
  
  # Computing the expected sufficient marginals 
  c_hat_k_1 <- function(beta,f){
    n <- length(n_r)
    n_beta <- nrow(beta)
    vec <- matrix(rep(0,n_beta*m),nrow=n_beta,ncol=m)

    for (x in 1:m){
      for (j in 1:n_beta){
        for (i in 1:n){
          vec[j,x] <- vec[j,x]+f(beta,n_r[i],r[i],j,x)
          
        }
      }}
      c(vec)
    }
    
  
  
  # Running algorithm 
  beta <- matrix(rep(1,m*n_item),ncol=m,nrow=n_item) #Initial value of betas are set to 1 
  C <- c_hat_k_1(beta,e_ri_1)
  mat_fun <- function(beta){
    mat <- matrix(ncol=n_item*m, nrow=n_score)
   for(x in 1:m){
     k <-ifelse(x==1,0,n_item*(x-1))
     for ( i in 1:n_item){
       for ( j in 1:n_score) {
         mat[j,k+i] <- e_ri_1(beta,n_r[j],r[j],i,x)
       }
       
     } 
   } 
    mat <- cbind(mat,rowSums(mat))
    mat <- rbind(mat,colSums(mat))
    mat
  }
  for (t in 1:maxiter){
    if (is.na(sum((C - W)^2) < epsilon))
    {warning("Model failure")
      beta <- matrix(rep(NA,m*n_item),ncol=m,nrow=n_item)
      break    
      }
    if(sum((C - W)^2) <= epsilon ){
      
      break
    }

      beta <- update_beta(W/C,beta)
      C <- c_hat_k_1(beta,e_ri_1)
    }
    
  
  structure(list("beta" = log(beta),
       "r" = r,
       "n_items" = n_item, 
       "m" = m,
       "V" = V,
       "x" = x,
       "W" = W,
       "data" =data,
       "i"= t, 
       "eps" = epsilon, 
       "table" = if(return_tab){mat_fun(beta)}else{NA}),
       class = "RACH_PCM_IPF")
  
}


print.RACH_PCM_IPF <- function(list){
  mat <- list$table
  dimnames(mat) <- list(c(list$r,"C"),c(names(list$W),"R"))
  dimnames(list$beta) <- list(1:list$n_items,1:list$m)
  cat(" \n",
      "Iterative proportional fitting polytomous Rasch model:\n Estimates and table \n
      Item difficulty parameters\n" )
  print(list$beta)
  cat("\n\n Levels: ", 1:list$m,
      "\n Items: ", 1:list$n_items)
  cat("\n\nTable:\n\n")
  print(mat)
  cat("\n\n Row target: ",list$V,
      "\n Column target: ", list$W)
  cat("\n \n Iterations:" ,list$i)
}




## Function for computing bootstraped confidence intervals and std. error. 

Std_error_BS_PCM <- function(obj, data, B=1000,n=NA){
  n<- ifelse(is.na(n),nrow(obj$data),n) # setting n allows for more draws than the sample size. 
  n_item <- ncol(data)
  m <- obj$m
  beta_mat <- matrix(ncol=n_item*m,nrow=B)
  for (i in 1:B){
    temp <- data[sample(nrow(data), replace = TRUE), ] #sample with replacement
    est <- CML_IPF_PCM(temp, removeExtreme = TRUE, epsilon = obj$eps, return_tab = FALSE)
    beta_mat[i,] <- est$beta
    print(i)
  }
  beta_mat <- beta_mat
  sd <- apply(beta_mat, 2, function(x){sd(x,na.rm=TRUE)})
  structure(list("beta_sim" = beta_mat,
                 "beta" = c(obj$beta),
                 "sd" = sd,
                 "B" = B), class = "IPF_rasch_CI_PCM")
}




print.IPF_rasch_CI_PCM <- function(list){
  cat(" \n",
      "Iterative proportional fitting polytomous Rasch model:\n Estimates and confidence intervals with method bootstrap \n") 
  print(matrix(c(list$beta,
                 list$sd,
                 list$beta-1.96*list$sd,
                 list$beta+1.96*list$sd),ncol=4,
               dimnames=list(rep(1:(length(list$beta)/2),2),
                             c("Estimate", "Std. Error", "lower CI", "upper CI"))))
  
  cat("\n Repititions:" ,list$B)
}



autoplot.IPF_rasch_CI_PCM <- function(obj, grid=FALSE) {
  n <- length(obj$beta)
  d_plots <- list()
  for (i in 1:n){
    d_plots[[i]] <- ggplot() + geom_histogram(aes(x=obj$beta_sim[,i],y=..density..), col="black", fill="orange") + 
      geom_vline(aes(xintercept= obj$beta[i]),color="red") + 
      geom_vline(aes(xintercept= mean(obj$beta_sim[,i],na.rm=TRUE))) +
      geom_vline(aes(xintercept=obj$beta-1.96*obj$sd), linetype="dashed")+ 
      geom_vline(aes(xintercept=obj$beta+1.96*obj$sd), linetype="dashed") + theme_bw()
    
  }
  
  if(grid){do.call("grid.arrange", c(d_plots, ncol=2))}
  
  
  else{
    for ( i in 1:n){
      print(d_plots[[i]])
    }}
}

