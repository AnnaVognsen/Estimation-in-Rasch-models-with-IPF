
CML_IPF_DIF <- function(data,items_dif,group,parinit=NA, 
                        maxiter=1000, 
                        epsilon=0.0001,
                        return_tab = TRUE) {
  data <- as.matrix(data)
  group <- as.vector(group)

  # Remove Na groups and extreme total scores
  if (any(is.na(group))){warning("NA groups and observations removed")}
  ind <- !(rowSums(data)==0|rowSums(data)==max(rowSums(data))|is.na(group))
  data <- data[ind,]
  Levels <- names(table(group))
  X <- 0:(length(Levels)-1)
  group <- as.integer(factor(group[ind]))-1

  # Data manipulation to get sufficient sums
  data_temp <- data[,-items_dif]
  nam <- (1:ncol(data))[-items_dif]
  colnames <- dimnames(data)[[2]]
  X_indicator<- rep(NA,ncol(data_temp))
  

  for(i in items_dif){
    for (x in X){
      data_temp <- cbind(data_temp,ifelse(group==x,data[,i],NA))
      nam <- c(nam,i)
      X_indicator <- c(X_indicator,x)
      
    }
  }

  # Par_indicator makes sure each group have seperate item difficulties
  par_indicator <- matrix(ncol=length(X), nrow=ncol(data_temp))
  for(i in 1:ncol(data_temp)){
    for (x in X){
      par_indicator[i,x+1] <- !(any(nam[i] == items_dif) & X_indicator[i]!=x)
    }
  }
  
  r_sums<- rowSums(data)
  colnames(data_temp) <- 1:ncol(data_temp)
  tab1 <- cbind(data_temp, Rowtotal = rowSums(data_temp, na.rm=TRUE))
  
  # Aggregate by the Rowtotal column
  grouped_data1 <- aggregate(tab1,FUN=function(x){sum(x,na.rm=TRUE)}, by=list(rowSums(data_temp, na.rm=TRUE)))
  
  
  dimnames(data_temp)[[2]]<- nam
  n_r_simp <- table(r_sums)
  r_sums_x <- cbind(r_sums,group)
  r <- as.integer(names(n_r_simp)) # Score in each scoregroup 
  n_r <- table(r_sums_x[,1],r_sums_x[,2])
  n_item <- ncol(data_temp)

  # Observed sufficient sums 
  W <- colSums(data_temp,na.rm=TRUE)          
  V <- n_r_simp*r 
  
  n_score <- length(r) 

  # Restrict item difficulties such that they sum to 0 
  
  prod1_DIF <- function(par) {
    par0 <- par[!duplicated(nam)]
    if(prod(par0)==1){par}
    else {
      geometric_mean <- exp(mean(log(par0)))
      par <- par / geometric_mean
      par
      
    }
  }

  # gamma function
  sum_prod <- function(r,delta){
    temp <- combn(delta,r,FUN=prod)
    sum(temp)
  }
  
  gamma_fun <- function(r_p,delta,k=ncol(data)){
    if(r_p==0){x<- 1}
    else if (r_p==1){x<- sum(delta)}
    else if(r_p==k){x<-1}
    else {x <- sum_prod(r_p,delta)}
    x
  }
  
  update_beta <- function(beta,X){
    beta_upt <- X*beta
    prod1_DIF(beta_upt)
    
  }


  # Expected table entries
  e_ri_1 <- function(delta,r,i) {
    col <- as.integer(rownames(n_r)) == r
    if(nam[i] %in% items_dif){
      x <- X_indicator[i]
      
      par <- delta[par_indicator[, x+1]]
      id <- names(par)==nam[i]
      n_r[col,x+1]*par[id]*(gamma_fun(r-1,par[!id])/gamma_fun(r,par))
    
      
      }
    else{
      temp <- 0
      for (x in X){
        
        par <- delta[par_indicator[, x+1]]
        temp<- temp + n_r[col,x+1]*par[i]*(gamma_fun(r-1,par[-i])/gamma_fun(r,par))}
        temp
        }
    
    }
    
  

# Expected sufficient sums
  c_hat_k_1 <- function(beta,f){
    n <- length(r)
    n_beta <- n_item
    vec <- numeric(n_beta)

    for (j in 1:n_beta){
        for (i in 1:n){
          vec[j] <- vec[j]+e_ri_1(beta,r[i],j)
        }
        
      }
    vec
  }
  
  # Running algorithm 
  delta <- rep(1,n_item)
  names(delta) <- nam
  C <- c_hat_k_1(delta,e_ri_1)

  #Expected contingency table 
  mat_fun <- function(beta){
    mat <- matrix(ncol=n_item, nrow=n_score)
    n <- length(r)
    n_beta <- length(beta)
    
   
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
      delta <- rep(NA,n_item)
      break    }
    if(sum((C - W)^2) <= epsilon * (sum(W^2) + epsilon)){
      
      break
    }
    
    delta <- update_beta(delta,W/C)
    C <- c_hat_k_1(delta,e_ri_1)
  }
 # getting estimated parameters
  par <- delta[!duplicated(nam)]
  par_star <- delta[!(delta%in% par)]
  beta <- par[names(par)==items_dif]/par_star
  structure(
    list("delta" = -log(par[sort(unique(nam),index.return=TRUE)$ix]),
         "beta" = log(beta),
         "nam" = nam,
       "r" = r,
       "Levels" = Levels, 
       "item_names" = colnames,
       "itemDIF" = items_dif,
       "V" = V,
       "W" = W,
       "data" =data,
       "group"=group,
       "i"= t, 
       "eps" = epsilon, 
       "mat_fun" = mat_fun,
       "data_tab" = grouped_data1[,c(-1,-(ncol(grouped_data1)))],
       "table" = if(return_tab){mat_fun(delta)} else{NA},
       "method" = "CML"), class="DIF_RM")
  
}


# Generic print function
print.DIF_RM <- function(list){
  
  mat <- list$table
  dimnames(mat) <- list(c(list$r,"C"),c(list$nam,"R"))
  
  n_item <- length(list$delta)
  items <- 1:n_item
  
  item_dif <- list$itemDIF
  X <- list$Levels
  no_dif <- rep(NaN, n_item)
  nam <- c("alpha")
  tab <-cbind(list$delta)
  for (i in item_dif){
    for(x in 1:length(X[-1])){
      temp <- no_dif
      temp[i] <- (list$beta[names(list$beta)==i])[x]
      tab <- cbind(tab, temp)
      nam <- cbind(nam,paste("beta_",X[1+x])) 
    }
  }
  dimnames(tab) <- list(list$item_names,nam)
  cat("\n",
      "Iterative proportional fitting Rasch model with DIF:\n Estimates and table \n Item difficulty parameters and DIF parameters \n" )
  print(tab)
  cat("Levels:",list$Levels, "\nIntercept level is ",list$Levels[1])
  cat("\nTable:\n")
  print(mat)
  cat("\n Row target: ",list$V,
      "\n Column target: ", list$W)
  cat("\n \n Iterations:" ,list$i,"\n Method:", list$method)
}




## Computing standard error 


Std_error_DIF<- function(obj, data, X, B=1000){
  delta <- obj$delta
  beta <- obj$beta
  
  item <- obj$itemDIF

  delta_mat <- matrix(ncol=length(delta),nrow=B)
  beta_mat <- matrix(ncol=length(beta), nrow=B)
  for (i in 1:B){
    id <- sample(nrow(data), replace = TRUE)
    temp <- data[id,] #sim_data(beta)
    x_temp <- X[id]
    est <- suppressWarnings(CML_IPF_DIF(temp,item,x_temp,epsilon = obj$eps, return_tab = FALSE))
    beta_mat[i,] <- tryCatch(est$beta, error = function(e){NA})
    delta_mat[i,] <- est$delta
    print(i)
  }

  sd_delta <- apply(delta_mat, 2, function(x){sd(x,na.rm=TRUE)})
  sd_beta <- apply(beta_mat, 2, function(x){sd(x,na.rm=TRUE)})
  structure(list("beta_sim" = beta_mat,
                 "delta_sim" = delta_mat,
                 "delta" = delta,
                 "beta" = beta,
                 "sd_beta" = sd_beta,
                 "sd_delta" = sd_delta,
                 "item_names" = obj$item_names,
                 "itemDIF" = obj$itemDIF,
                 "B" = B), class = "IPF_rasch_CI_DIF")
}

#generic print function for standard errors 

print.IPF_rasch_CI_DIF <- function(list){
  cat(" \n",
      "Iterative proportional fitting dichotomous Rasch model with DIF:\n Estimates and confidence intervals with method bootstrap \n") 
  print(matrix(c(list$delta, 
                 list$beta,
                 list$sd_delta,
                 list$sd_beta,
                 list$delta-1.96*list$sd_delta,
                 list$beta-1.96*list$sd_beta,
                 list$delta+1.96*list$sd_delta,
                 list$beta+1.96*list$sd_beta),ncol=4,
               dimnames=list(c(paste("alpha_",list$item_names),paste("Beta_",list$item_names[list$itemDIF])),
                             c("Estimate", "Std. Error", "lower CI", "upper CI"))))
  
  cat("\n Repititions:" ,list$B)
}

# autoplot function 

autoplot.IPF_rasch_CI_DIF<- function(obj, grid=FALSE) {
  n <- length(obj$delta) 
  n_beta <- length(obj$beta)
  d_plots <- list()
  
  delta_lwr <- obj$delta-1.96*obj$sd_delta
  beta_lwr <- obj$beta-1.96*obj$sd_beta
  delta_upr <- obj$delta+1.96*obj$sd_delta
  beta_upr <- obj$beta+1.96*obj$sd_beta
  
  for (i in 1:n){
    d_plots[[i]] <- ggplot() + geom_histogram(aes(x=obj$delta_sim[,i],y=..density..), col="black", fill="orange") + 
      geom_vline(aes(xintercept= obj$delta[i]),color="red") + 
      geom_vline(aes(xintercept= mean(obj$delta_sim[,i],na.rm=TRUE))) +
      geom_vline(aes(xintercept=delta_lwr[i]), linetype="dashed")+ 
      geom_vline(aes(xintercept=delta_upr[i]), linetype="dashed") + theme_bw()
    
  }
  if (ncol(obj$beta_sim) > 1) {
    for(j in (n+1):(n+n_beta)){
      i <- j - n
      d_plots[[j]] <- ggplot() + geom_histogram(aes(x=obj$beta_sim[,i],y=..density..), col="black", fill="orange") + 
        geom_vline(aes(xintercept= obj$beta[i]),color="red") + 
        geom_vline(aes(xintercept= mean(obj$beta_sim[,i],na.rm=TRUE))) +
        geom_vline(aes(xintercept=beta_lwr[i]), linetype="dashed")+ 
        geom_vline(aes(xintercept=beta_upr[i]), linetype="dashed") + theme_bw()
    }
  } else {
    d_plots[[n+1]] <- ggplot() + geom_histogram(aes(x=obj$beta_sim,y=..density..), col="black", fill="orange") + 
      geom_vline(aes(xintercept= obj$beta),color="red") + 
      geom_vline(aes(xintercept= mean(obj$beta_sim,na.rm=TRUE))) +
      geom_vline(aes(xintercept=beta_lwr), linetype="dashed")+ 
      geom_vline(aes(xintercept=beta_upr), linetype="dashed") + theme_bw()
  }
  
  
  if(grid){do.call("grid.arrange", c(d_plots, ncol=2))}

  else{
    for (i in 1:(n+n_beta)){
      print(d_plots[[i]])
    }}
  return(d_plots)
}



