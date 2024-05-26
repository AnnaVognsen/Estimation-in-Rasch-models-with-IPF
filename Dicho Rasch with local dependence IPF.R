


  
  
CML_IPF_LD <- function(data,items_LD,parinit=NA, 
                         maxiter=1000, 
                         epsilon=0.0001,
                         return_tab = TRUE) {
    if (length(items_LD)>2) {stop("only 2 locally dependent items aloud")}
    data <- as.matrix(data)
    ind <- !(rowSums(data)==0|rowSums(data)==ncol(data))) # removing extreme total scores
    data <- data[ind,]
    interceptLD <- items_LD[1] # first dependent item is intercept 
    Levels <- names(table(data[,interceptLD]))
    X <- 0:(length(Levels)-1) # values first item can take
    group <- data[,interceptLD]
    i_dep <- items_LD[2] 


  #Data manipulations to get sufficient sums 
    data_temp <- data[,-i_dep]
    nam <- (1:ncol(data))[-i_dep]
    colnames <- dimnames(data)[[2]]
    X_indicator<- rep(NA,ncol(data_temp))

    
    for(i in i_dep){
      for (x in X){
        data_temp <- cbind(data_temp,ifelse(group==x,data[,i],NA))
        nam <- c(nam,i)
        X_indicator <- c(X_indicator,x)
        
      }
    }

  # Parameter indicator for making sure those with y_i=1 fase one set of item difficulties and y_i=0 fase others. 
  
    par_indicator <- matrix(ncol=length(X), nrow=ncol(data_temp))
    for(i in 1:ncol(data_temp)){
      for (x in X){
        par_indicator[i,x+1] <- !(any(nam[i] == i_dep) & X_indicator[i]!=x)
      }
    }

  # Observed contigency table 
    tab1 <- cbind(data_temp, Rowtotal = rowSums(data_temp, na.rm=TRUE))
    # Aggregate by the Rowtotal column
    grouped_data1 <- aggregate(tab1,FUN=function(x){sum(x,na.rm=TRUE)}, by=list(rowSums(data_temp, na.rm=TRUE)))

    dimnames(data_temp)[[2]]<- nam
    r_sums<- rowSums(data)
    n_r <- table(r_sums)
    r <- as.integer(names(n_r)) # Score in each scoregroup 
    n_item <- ncol(data_temp)

)
    
  # target table margins 
    W <- colSums(data_temp,na.rm=TRUE)          
    V <- n_r*r 
    
  # number of rows 
    n_score <- length(r) 

  # Restriction of item difficulties 
    prod1<- function(par) {
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
      else if(r_p==k){x<-prod(delta)}
      else if(r_p<0|length(delta)==0){x <- 0}
      else {x <- sum_prod(r_p,delta)}
      x
    }
    
    gamma_fun_LD <- function(r_p,delta0,delta1,k=ncol(data)){
      if(r_p==0){x<- 1}
      else if (r_p==1){x<- sum(delta0)}
      else if(r_p==k){x<-prod(delta1)}
      
      
      else {
        id <- names(delta1)==interceptLD
        x <- delta1[id]*gamma_fun(r_p-1,delta1[!id])+gamma_fun(r_p,delta0[!id])
      
      }
      x
    }
    
    update_beta <- function(beta,X){
      beta_upt <- X*beta
      prod1(beta_upt)
      
    }

  #expected entries 
    e_ri_1 <- function(delta,r,i) {
      
      if(nam[i] %in% i_dep){
        x <- X_indicator[i]
        par <- delta[par_indicator[, x+1]]
        par1 <- delta[par_indicator[, 2]]
        par0 <- delta[par_indicator[, 1]]
        id <- names(par)==nam[i]
        par_inter <- ifelse(x==0,1,par[names(par)==interceptLD])
       
        n_r[r]*par_inter*par[id]*(gamma_fun(r-1-x,par[!id & names(par)!=interceptLD])/gamma_fun_LD(r,par0,par1))
        
        
      }
      
      else{
        if(nam[i] %in% interceptLD){
          par1 <- delta[par_indicator[, 2]]
          par0 <- delta[par_indicator[, 1]]
          
          n_r[r]*par1[i]*(gamma_fun(r-1,par1[-i])/gamma_fun_LD(r,par0,par1))

        }
        else{
          par1 <- delta[par_indicator[, 2]]
          par0 <- delta[par_indicator[, 1]]
          n_r[r]*par0[i]*(gamma_fun_LD(r-1,par0[-i],par1[-i])/gamma_fun_LD(r,par0,par1))
          
          
        } 
      }
      
      
    }
    
    
    
    # Expected column sums
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
      if (is.na(sum((C - W)^2) < epsilon ))
      {warning("Model failure")
        delta <- rep(NA,n_item)
        break    }
      if(sum((C - W)^2) <= epsilon * (sum(W^2) + epsilon)){
        
        break
      }
      
      delta <- update_beta(delta,W/C)
      C <- c_hat_k_1(delta,e_ri_1)
    }

  # getting parameters 
    par <- delta[!duplicated(nam)]
    par_star <- delta[!(delta%in% par)]
    beta <- par[names(par)==i_dep]/par_star
    structure(
      list("delta" = -log(par[sort(unique(nam),index.return=TRUE)$ix]),
           "zeta"= -log(beta),
           "data_tab"=grouped_data1,
           "r" = r,
           "Levels" = Levels, 
           "item_names" = colnames,
           "itemDEP" = items_LD,
           "mat_fun"=mat_fun,
           "nam" = nam,
           "V" = V,
           "W" = W,
           "data" =data,
           "i"= t, 
           "eps" = epsilon,
           "method" = "CML",
           "table" = if(return_tab){mat_fun(delta)} else{NA}), class="LD_RM")
    
}

# A generic print function
print.LD_RM <- function(list){
  mat <- if(is.matrix(list$table)){list$table} else{"Not calculated"}
  dimnames(mat) <- list(c(list$r,"C"),c(names(list$delta),names(list$zeta),"R"))
  
  n_item <- length(list$delta)
  items <- 1:n_item
  
  item_dep <- list$itemDEP

  nam <- c("alpha")
  tab <-list$delta
  

  cat("\n",
      "Iterative proportional fitting Rasch model with local dependence: \n Item difficulty parameters  \n" )
  print(tab)
  cat("\nDependence parameter beteen items ",item_dep, "\nzeta ",list$zeta)
  
  cat("\n\nTable:\n")
  print(mat)
  cat("\n Row target: ",list$V,
      "\n Column target: ", list$W)
  cat("\n \n Iterations:" ,list$i,"\n Method:", list$method)
}





#Bootstrapping standard error and confidence intervals for parameters

Std_error_LD<- function(obj, data, X, B=1000){
  delta <- obj$delta
  zeta <- obj$zeta
  item <- obj$itemDEP
  
  delta_mat <- matrix(ncol=length(delta),nrow=B)
  zeta_mat <- matrix(ncol=length(zeta), nrow=B)
  for (i in 1:B){
    id <- sample(nrow(data), replace = TRUE)
    temp <- data[id,] #sim_data(beta)
    est <- suppressWarnings(CML_IPF_LD(temp,item,epsilon = obj$eps, return_tab = FALSE))
    zeta_mat[i,] <- ifelse(length(est$zeta)==0,NA,est$zeta)
    delta_mat[i,] <- est$delta
    print(i)
  }
  
  sd_delta <- apply(delta_mat, 2, function(x){sd(x,na.rm=TRUE)})
  sd_zeta <- apply(zeta_mat, 2, function(x){sd(x,na.rm=TRUE)})
  structure(list("zeta_sim" = zeta_mat,
                 "delta_sim" = delta_mat,
                 "delta" = delta,
                 "zeta" = zeta,
                 "sd_zeta" = sd_zeta,
                 "sd_delta" = sd_delta,
                 "item_names" = obj$item_names,
                 "itemDEP" = obj$itemDEP,
                 "B" = B), class = "IPF_rasch_CI_LD")
}

# Generic print functions
print.IPF_rasch_CI_LD <- function(list){
  cat(" \n",
      "Iterative proportional fitting dichotomous Rasch model with local dependence:\n Estimates and confidence intervals with method bootstrap \n") 
  print(matrix(c(list$delta, 
                 list$zeta,
                 list$sd_delta,
                 list$sd_zeta,
                 list$delta-1.96*list$sd_delta,
                 list$zeta-1.96*list$sd_zeta,
                 list$delta+1.96*list$sd_delta,
                 list$zeta+1.96*list$sd_zeta),ncol=4,
               dimnames=list(c(paste("delta_",list$item_names),paste("zeta_",sep="",list$itemDEP[1],",",list$itemDEP[2])),
                             c("Estimate", "Std. Error", "lower CI", "upper CI"))))
  
  cat("\n Repititions:" ,list$B)
}


# Autoplot
autoplot.IPF_rasch_CI_LD<- function(obj, grid=FALSE) {
  n <- length(obj$delta) 
  n_zeta <- length(obj$zeta)
  d_plots <- list()
  delta_lwr <- obj$delta-1.96*obj$sd_delta
  zeta_lwr <- obj$zeta-1.96*obj$sd_zeta
  delta_upr <- obj$delta+1.96*obj$sd_delta
  zeta_upr <- obj$zeta+1.96*obj$sd_zeta
  for (i in 1:n){
    d_plots[[i]] <- ggplot() + geom_histogram(aes(x=obj$delta_sim[,i],y=..density..), col="black", fill="orange") + 
      geom_vline(aes(xintercept= obj$delta[i]),color="red") + 
      geom_vline(aes(xintercept= mean(obj$delta_sim[,i],na.rm=TRUE))) +
      geom_vline(aes(xintercept=delta_lwr), linetype="dashed")+ 
      geom_vline(aes(xintercept=delta_upr), linetype="dashed") + theme_bw()
    
  }
  if (ncol(obj$zeta_sim) > 1) {
    for(j in (n+1):(n+n_zeta)){
      i <- j - n
      d_plots[[j]] <- ggplot() + geom_histogram(aes(x=obj$zeta_sim[,i],y=..density..), col="black", fill="orange") + 
        geom_vline(aes(xintercept= obj$zeta[i]),color="red") + 
        geom_vline(aes(xintercept= mean(obj$zeta_sim[,i],na.rm=TRUE))) +
        geom_vline(aes(xintercept=zeta_lwr), linetype="dashed")+ 
        geom_vline(aes(xintercept=zeta_upr), linetype="dashed") + theme_bw()
    }
  } else {
    d_plots[[n+1]] <- ggplot() + geom_histogram(aes(x=obj$zeta_sim,y=..density..), col="black", fill="orange") + 
      geom_vline(aes(xintercept= obj$zeta),color="red") + 
      geom_vline(aes(xintercept= mean(obj$beta_sim,na.rm=TRUE))) +
      geom_vline(aes(xintercept=zeta_lwr), linetype="dashed")+ 
      geom_vline(aes(xintercept=zeta_upr), linetype="dashed") + theme_bw()
  }
  
  
  if(grid){do.call("grid.arrange", c(d_plots, ncol=2))}
  
  else{
    for (i in 1:(n+n_beta)){
      print(d_plots[[i]])
    }}
}
