# function for calculating likelihood ratio test and p-values 
# obj1 is the null model 
# obj2 is the alternative model 
anova_IPF <- function(obj1,obj2){
  if(class(obj1)== "Rasch_IPF_Estimation" & class(obj2) == "DIF_RM"){
    # Anova for DIF in dichotomous rasch model 
    
    delta_1 <- exp(-obj1$delta[obj2$nam])
    
    names(delta_1) <- obj2$nam
    tab1 <- obj2$mat_fun(delta_1)
    tab1 <- tab1[-nrow(tab1),-ncol(tab1)]
    tab2 <- obj2$table
    tab2 <- tab2[-nrow(tab2),-ncol(tab2)]
    
    LRT <- 2*sum(tab1*log(tab1/tab2),na.rm=TRUE)
    df <- length(obj2$beta)
    p <- pchisq(LRT,df,lower.tail = FALSE)
    
  }
  
  else if (class(obj1)== "Rasch_IPF_Estimation" & class(obj2) == "LD_RM"){
    # Anova for dichotomous rasch model with local dependence 
    delta_1 <- exp(-obj1$delta[obj2$nam])
    names(delta_1) <- obj2$nam
    tab1 <- obj2$mat_fun(delta_1)
    tab1 <- tab1[-nrow(tab1),-ncol(tab1)]
    tab2 <- obj2$table
    tab2 <- tab2[-nrow(tab2),-ncol(tab2)]
    
    LRT <- 2*sum(tab1*log(tab1/tab2),na.rm=TRUE)
    df <- length(obj2$zeta)
    p <- pchisq(LRT,df,lower.tail = FALSE)
    
    
  }
  
  
  else if(class(obj1)== "RACH_PCM_IPF" & class(obj2) == "RACH_PCM_IPF_DIF") {
    # Anova for polytomous Rasch with DIF 
    delta1 <- c()
    for(i in 1:length(obj2$X_indicator)){
      
        delta1 <- c(delta1,obj1$beta[obj2$nam[i],obj2$X_indicator[i]])
      
    }
    delta1 <- exp(delta1)
    names(delta1) <- obj2$nam
    tab1 <- obj2$mat_fun(delta1)
    tab1 <- tab1[-nrow(tab1),-ncol(tab1)]
    tab2 <- obj2$table
    tab2 <- tab2[-nrow(tab2),-ncol(tab2)]
    
    LRT <- 2*sum(tab1*log(tab1/tab2),na.rm=TRUE)
    
    df <- length(obj2$beta_x)
    p <- pchisq(LRT,df,lower.tail = FALSE) 
    
  }
  
  else{
    stop("obj1 must be class Rasch_IPF_Estimation or RACH_PCM_IPF and obj2 must be DIF_RM or RACH_PCM_IPF_DIF")
  }
  
  list("LRT"=LRT,"p"=p,"df"=df)#cat("Null hypothesis: beta_x=0","P-value:",p, "LRT:",LRT, "degrees of freedom: ",df )
}

