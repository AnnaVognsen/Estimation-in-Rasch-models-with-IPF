CML_IPF_PCM_DIF <- function(data,
                            items_dif,
                            group,
                            removeExtreme = TRUE,
                            x=NA,
                            parinit=NA, 
                            maxiter=1000, 
                            epsilon=0.0001 ,
                            return_tab = TRUE) {
  data <- as.matrix(data)
  group <- as.vector(group)
  if(is.na(x)){
    m <- max(data)
    x <- seq(1,m)}
  
  Levels <- names(table(group)) # get the levels from the group variable before it is transformed into 0,1,2
  n_item <- ncol(data) # number of columns = items
  
  if(removeExtreme){ # remove itemsresponses with r=0 or r=m*k also from data and group variable. 
    ind <- !(rowSums(data)==0|rowSums(data)==m*n_item|is.na(group))
    data <- data[ind,]
    group <- group[ind]
  } 
  group <- as.integer(factor(group))-1 # group to integer
  G <- 0:(length(Levels)-1)
  
  
  ## Data manipulations for getting the observed sufficient sums 
  
  data_temp_1 <- (data==x[1])*x[1]
  colnames <- rep(dimnames(data)[[2]],m)
  nam <- rep(1:n_item,m)
  for (xhat in x[-1]){
    data_temp_1 <- cbind(data_temp_1,(data==xhat)*xhat)
    
  }
  
  data_temp <- data_temp_1[,nam!=items_dif]
  colnames(data_temp_1) <- nam
  nam <- nam[nam!=items_dif]
  
  
  ## Creating indicator vectors with an entree for each collumn for each sufficient sum. 
  # X Indicates y_i = 1 or y_i = 2
  # G indicates which groups are considered in the collumn 
  # Nam indicates which item sum we consider
  # Par_indicator has a column for each group in G. Each group have different item difficulties. Par_indicator keeps track of this 
  X_indicator <- rep(1,ncol(data_temp)/2)
  for (xhat in x[-1]){
    X_indicator <- c(X_indicator,rep(xhat,ncol(data_temp)/2))
  }  
  
  G_indicator<- rep(NA,ncol(data_temp))
  
  
  for(i in items_dif){
    for (g in G){
      for(xhat in x){
        temp <- data_temp_1[,colnames(data_temp_1) == i]
        data_temp <- cbind(data_temp,ifelse(group==g,temp[,xhat],NA))
        nam <- c(nam,i)
        G_indicator <- c(G_indicator,g)
        X_indicator <- c(X_indicator,xhat)
      }
      
      
    }
  }
  
  par_indicator <- matrix(ncol=length(G), nrow=ncol(data_temp))
  for(i in 1:ncol(data_temp)){
    for (g in G){
      par_indicator[i,g+1] <- !(any(nam[i] == items_dif) & G_indicator[i]!=g)
    }
  }
  
  r_sums <- rowSums(data)
  n_r_simple <- table(r_sums)
  r <- as.integer(names(n_r_simple)) # Score in each scoregroup 
  n_r_simple <- unname(n_r_simple)
  n_score <- length(r) 
  r_sums_x <- cbind(r_sums,group)
  n_r <- table(r_sums_x[,1],r_sums_x[,2]) # count persons with R=r and G=g
  colnames(data_temp) <-nam
 
  
  # Observed sufficient margins 
  W <- colSums(data_temp,na.rm=TRUE)
  V <- n_r_simple*r
  
  
  ## calculate the observed data table for anova test
  tab1 <- data.frame(cbind(data_temp, Rowtotal = rowSums(data_temp, na.rm=TRUE)))
  # Aggregate by the Rowtotal column
  grouped_data1 <- aggregate(tab1,FUN=function(x){sum(x,na.rm=TRUE)}, by=list(tab1$Rowtotal))
  
  
  # Gamma recurisive funciton 
  gamma_rec <- memoise::memoise(function(beta,r){
    beta <- matrix(beta,ncol=m)
    k<- n_item
    if (r==0){return(1)}
    if (r==1){return(sum(beta[,1]))}
    if (r==(k*m)){return(prod(beta[,m]))}
    if(r<0|nrow(beta)==0|r>(nrow(beta)*k)){return(0)}
    
    else{
      i <- 1
      temp <- 0 
      for (xhat in 1:m){
        temp <- temp + beta[i,xhat]*gamma_rec(beta[-i,],r-xhat)
      }
      return(temp + gamma_rec(beta[-i,],r))        
    }
    
  })
  

  ## prod1_pc makes sure the item difficulties sum to 0. 
  prod1_pc <- function(matrix) {
    par0 <- matrix[is.na(G_indicator)|G_indicator==0]
    if(prod(par0)==1){return(matrix)}
    else {
      geometric_mean <- exp(mean(log(par0)))
      matrix <- (matrix / geometric_mean)
      
      
      return(matrix)
    }
    
  }
  
  
  
  
  update_beta<-function(X,beta){
    beta_upt <- X*as.vector(beta)
    
    prod1_pc(beta_upt)
  }
  
  
  
  ## Expected cell frequencies. 
  e_ri_1 <- function(delta,r,i) {
    if(nam[i] %in% items_dif){
      g <- G_indicator[i]
      xhat <- X_indicator[i]
      p <- delta[par_indicator[, g+1]]
      par <- rbind(matrix(p[1:((n_item-length(items_dif))*m)],ncol=m),
                   p[-(1:((n_item-length(items_dif))*m))])
      id <- nam[i]
      n_r[r,g+1]*par[id,xhat]*xhat*(gamma_rec(par[-id,],r-xhat)/gamma_rec(par,r))
      
      
    }
    else{
      temp <- 0
      xhat <- X_indicator[i]
      for (g in G){
        id <- nam[i]
        p <- delta[par_indicator[, g+1]]
        par <- rbind(matrix(p[1:((n_item-length(items_dif))*m)],ncol=m),p[-(1:((n_item-length(items_dif))*m))])
        temp<- temp + n_r[r,g+1]*par[id,xhat]*xhat*(gamma_rec(par[-id,],r-xhat)/gamma_rec(par,r))}
      temp
    }
    
  }
  
  ## Calculates the expected sufficient marginals. 
  c_hat_k_1 <- function(beta,f){
    n <- length(r)
    n_beta <- length(beta)
    vec <- numeric(n_beta)
    
    
    for (j in 1:n_beta){
      for (i in 1:n){
        vec[j] <- vec[j]+f(beta,r[i],j)
        
      }
    }
    c(vec)
  }
  

  # Running algorithm 
  i <- 1
  beta <- rep(1,ncol(data_temp))
  names(beta) <- nam
  C <- c_hat_k_1(beta,e_ri_1)
  mat_fun <- function(beta){
    n <- length(r)
    n_beta <- length(beta)
    mat <- matrix(ncol=n_beta*r, nrow=n)
    for (j in 1:n_beta){
      for (i in 1:n){
        mat[i,j] <- e_ri_1(beta,r[i],j)
        
      }
    }
    mat <- cbind(mat,rowSums(mat))
    mat <- rbind(mat,colSums(mat))
    mat
  }
  for (t in 1:maxiter){
    if (is.na(sum((C - W)^2) <= epsilon * (sum(W^2) + epsilon)))
    {warning("Model failure")
      beta <- rep(rep(NA,ncol(data_temp)),ncol=m)
      break    
    }
    if(sum((C - W)^2) <= epsilon ){
      
      break
    }
    
    
    beta <- update_beta(W/C,beta)
    C <- c_hat_k_1(beta,e_ri_1)
  }
  ## Transforming the resulting parameter vector into (alphaa_1,..,alpha_k) and (beta_X_1,..)
  p <- beta[par_indicator[,1]]
  p_star <- beta[!(beta%in%p)]
  beta_x <- p[names(p)==items_dif]/p_star
  beta_mat <- rbind(matrix(p[1:((n_item-length(items_dif))*m)],ncol=m),
                    p[-(1:((n_item-length(items_dif))*m))])
  
  
  
  structure(list("beta" = log(beta_mat),
                 "beta_x" = log(beta_x),
                 "ItemDIF"= items_dif,
                 "Levels" = Levels, 
                 "nam"= nam,
                 "r" = r,
                 "n_items" = n_item, 
                 "itemNames"=colnames,
                 "data_tab" = grouped_data1[,c(-1,-(ncol(grouped_data1)))],
                 "mat_fun" = mat_fun,
                 "m" = m,
                 "V" = V,
                 "x" = x,
                 "X_indicator" = X_indicator,
                 "W" = W,
                 "data" =data,
                 "i"= t, 
                 "eps" = epsilon, 
                 "table" = if(return_tab){mat_fun(beta)}else{NA},
                 "method" = "CML"),
            class = "RACH_PCM_IPF_DIF")
  
}



print.RACH_PCM_IPF_DIF <- function(list){
  
  mat <- list$table
  Test <- Test_IPF(list, method="Pearson")
  dimnames(mat) <- list(c(list$r,"C"),c(names(list$W),"R"))
  dimnames(list$beta) <- list(1:list$n_items,1:list$m)
  cat(" \n",
      "Iterative proportional fitting polytomous Rasch model with DIF:\n Estimates and table \n
      Item difficulty parameters\n" )
  print(list$beta)
  cat("\n\nOrdinal Levels: ", 1:list$m,
      "\nItems: ", 1:list$n_items)
  cat("\n\n\nItems with DIF:", list$ItemDIF,"\nGroup levels:",list$Levels,"\nIntercept",list$Levels[1] ,"\nBeta_x",list$beta_x )
  
  cat("\nGoodness of fit: Chi-squared test\n Method:",Test$method, "\nTest size:" ,Test$test_size, " p-value:", Test$p_value, "df:", Test$df)
  cat("\n\nTable:\n\n")
  print(mat)
  cat("\n\nRow target: ",list$V,
      "\nColumn target: ", list$W)
  cat("\n \nIterations:" ,list$i)
}



## Computing std error using bootstrap

Std_error_BS_PCM_DIF <- function(obj, data, X,B=200,n=NA){
  n<- ifelse(is.na(n),nrow(obj$data),n)
  n_item <- obj$n_items
  m <- obj$m
  beta_mat <- matrix(ncol=n_item*m,nrow=B)
  beta_x_mat <- matrix(ncol=length(obj$ItemDIF)*m,nrow=B)
  for (i in 1:B){
    print(i)
    samp <- sample(nrow(data),n, replace = TRUE)
    group <- X[samp]
    temp <- data[samp, ] #sim_data(beta)
    est <- tryCatch(CML_IPF_PCM_DIF(temp, obj$ItemDIF,
                                    as.vector(group),
                                    removeExtreme = TRUE, epsilon = obj$eps, return_tab = FALSE), 
                    error = function(A){NA}, 
                    warning = function(A){NA})
    if (is.list(est)){
      beta_mat[i,] <- est$beta
      beta_x_mat[i,] <- ifelse(est$beta_x==0,rep(NA,length(obj$ItemDIF)*m),est$beta_x) 
    }
  
  }
  
  sd_beta <- apply(beta_mat, 2, function(x){sd(x,na.rm=TRUE)})
  sd_beta_x <- apply(beta_x_mat, 2, function(x){sd(x,na.rm=TRUE)})
  structure(list("beta_sim" = beta_mat,
                 "beta_x_sim"=beta_x_mat,
                 "beta" = obj$beta,
                 "beta_x"=obj$beta_x,
                 "sd_beta" = sd_beta,
                 "sd_beta_x"=sd_beta_x,
                 "B" = B), class = "IPF_rasch_CI_PCM_DIF")
}

print.IPF_rasch_CI_PCM_DIF <- function(list){
  cat(" \n",
      "Iterative proportional fitting polytomous Rasch model with DIF:\n Estimates and confidence intervals with method bootstrap \n") 
  print(matrix(c(list$beta,
                 list$beta_x,
                 list$sd_beta,
                 list$sd_beta_x,
                 list$beta-1.96*list$sd_beta,
                 list$beta_x-1.96*list$sd_beta_x,
                 list$beta+1.96*list$sd_beta,
                 list$beta_x+1.96*list$sd_beta_x),ncol=4,
               dimnames=list(c(rep(1:nrow(list$beta),2),paste(names(list$beta_x),"_dif")),
                             c("Estimate", "Std. Error", "lower CI", "upper CI"))))
  
  cat("\n Repititions:" ,list$B)
}



