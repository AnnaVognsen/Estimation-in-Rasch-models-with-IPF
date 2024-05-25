CML_IPF_PCM_LD<- function(data,
                            items_LD,
                            removeExtreme = TRUE,
                            x=NA,
                            parinit=NA, 
                            maxiter=1000, 
                            epsilon=0.0001 ,
                            return_tab = TRUE) {
  
  if(is.na(x)){
    m <- max(data)
    x <- seq(1,m)}
       
  if(removeExtreme){
    ind <- !(rowSums(data)==0|rowSums(data)==max(rowSums(data)))
    data <- data[ind,]
  }
  
  i_dep <- items_LD[2]
  interceptLD <- items_LD[1]
  n_item <- ncol(data) # number of columns = items
  group <- data[,interceptLD]
  Levels <- names(table(group))
  G <- 0:(length(Levels)-1)
  
 
  
  data_temp_1 <- (data==x[1])*x[1]
  
  colnames <- rep(dimnames(data)[[2]],m)
  nam <- rep(1:n_item,m)
  for (xhat in x[-1]){
    data_temp_1 <- cbind(data_temp_1,(data==xhat)*xhat)
    
  }
  
  colnames(data_temp_1) <- nam
  data_temp <- data_temp_1[,nam!=i_dep]

  nam <- nam[nam!=i_dep]
  
  X_indicator <- rep(1,ncol(data_temp)/2)
  for (xhat in x[-1]){
    X_indicator <- c(X_indicator,rep(xhat,ncol(data_temp)/2))
  }  
  
  G_indicator<- rep(NA,ncol(data_temp))
  
  
  for(i in i_dep){
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
      par_indicator[i,g+1] <- !(any(nam[i] == i_dep) & G_indicator[i]!=g)
    }
  }
  

  index <- function(i){
    j<- length(x)
    if(i<j){return(i+1)}
    if(i==j){return(1)}
  }
  
  n <- nrow(data)
  r_sums <- rowSums(data)
  n_r <- table(r_sums)
  r <- as.integer(names(n_r)) # Score in each scoregroup 
  n_r <- unname(n_r)
  n_score <- length(r) 

  colnames(data_temp) <-nam
  # target table margins 
  W <- colSums(data_temp,na.rm=TRUE)
  V <- n_r*r 
  
  tab1 <- data.frame(cbind(data_temp, Rowtotal = rowSums(data_temp, na.rm=TRUE)))
  # Aggregate by the Rowtotal column
  grouped_data1 <- aggregate(tab1,FUN=function(x){sum(x,na.rm=TRUE)}, by=list(tab1$Rowtotal))
  
  gamma_rec <- memoise::memoise(function(beta,r){
    #browser()
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
  
  gamma_fun_LD <- function(r_p,delta0,delta1,delta2,k=n_item){
    if(r_p==0){x<- 1}
    else if (r_p==(k*m)){return(prod(delta2[,m]))}
    else if(r_p<0|nrow(delta0)==0|nrow(delta1)==0|nrow(delta2)==0){x<-0}
 
    
    else {
      #browser()
      id <- interceptLD
      
      x <- delta2[id,2]*gamma_rec(delta2[-id,],r_p-2)+delta1[id,1]*gamma_rec(delta1[-id,],r_p-1)+gamma_rec(delta0[-id,],r_p)
      
    }
    x
  }
  
  prod1_pc <- function(matrix) {
    par0 <- matrix[is.na(G_indicator)|G_indicator==0]
    if(prod(par0)==1){return(matrix)}
    else {
      geometric_mean <- exp(mean(log(par0)))
      matrix <- (matrix / geometric_mean)
      
      
      return(matrix)
    }
    
  }
  
  get_par <- function(p){
    rbind(matrix(p[1:((n_item-length(i_dep))*m)],ncol=m),
                  p[-(1:((n_item-length(i_dep))*m))])
  }
  
  
  update_beta<-function(X,beta){
    #browser()
    beta_upt <- X*as.vector(beta)
    
    prod1_pc(beta_upt)
  }
  
  
  e_ri_1 <- function(delta,r,i) {
    #browser()
    xhat <- X_indicator[i]
    
    p2 <- delta[par_indicator[, 3]]
    par2 <- get_par(p2)
    
    p1 <- delta[par_indicator[, 2]]
    par1 <- get_par(p1)
    
    p0 <- delta[par_indicator[, 1]]
    par0 <- get_par(p0)
    
    if(nam[i] %in% i_dep){
      #browser()
      g <- G_indicator[i]
      p <- delta[par_indicator[,g+1]]
      par <- get_par(p)
      
      
      id <- nam[i]
      par_inter <- ifelse(g==0,1,par[interceptLD,g])
      n_r[r]*par_inter*par[id,xhat]*xhat*(gamma_rec(par[-c(id,interceptLD),],r-xhat-g)/gamma_fun_LD(r,par0,par1,par2))
      #g?
  
      
    }
    
    else if(nam[i] %in% interceptLD){
      #browser()
      id <- nam[i]
      if(xhat==1){par_i <- par1} else{par_i <- par2}
      n_r[r]*par_i[id,xhat]*xhat*(gamma_rec(par_i[-id,],r-xhat)/gamma_fun_LD(r,par0,par1,par2))
      
    }
    
    
    else{
      id <- nam[i]
      n_r[r]*par0[id,xhat]*xhat*(gamma_fun_LD(r-xhat,par0[-id,],par1[-id,],par2[-id,])/gamma_fun_LD(r,par0,par1,par2))
      
      
    }
    
  }

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
  beta <- 1:ncol(data_temp)#rep(1,ncol(data_temp))
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
  #browser()
  for (t in 1:maxiter){
    if (is.na(sum((C - W)^2) <= epsilon * (sum(W^2) + epsilon)))
    {warning("Model failure")
      #beta <- rep(rep(NA,ncol(data_temp)),ncol=m)
      break    
    }
    if(sum((C - W)^2) <= epsilon ){
      
      break
    }
    
    
    beta <- update_beta(W/C,beta)
    C <- c_hat_k_1(beta,e_ri_1)
  }
  p <- beta[par_indicator[,1]]
  p_star <- beta[!(beta%in%p)]
  beta_x <- p[names(p)==i_dep]/p_star
  beta_mat <- rbind(matrix(p[1:((n_item-length(i_dep))*m)],ncol=m),
                    p[-(1:((n_item-length(i_dep))*m))])
  
  
  
  structure(list("beta" = log(beta_mat),
                 "zeta_x" = log(beta_x),
                 "ItemDIF"= i_dep,
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
            class = "RACH_PCM_IPF_LD")
  
}

