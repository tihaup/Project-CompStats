#-----------------------------------------------------------------------------
#### Data Generating Process
#-----------------------------------------------------------------------------

dgp <- function(p,N,treatment,regressors) { 
    
    if (regressors == "nocorr"){
        X <- matrix(rnorm(n=N*p,0,1), N, p)
    }
    
    if (regressors == "poly1"){
        X <- matrix(rnorm(n=N*p,0,1), N, p)
        X[,p-1] <- X[,1]^2
        X[,p] <- X[,2]^2
    }
    
    if (regressors == "poly2"){
        X <- matrix(rnorm(n=N*p,0,1), N, p)
        X[,2] <- X[,1]^2
    }    
    
    if (regressors == "interaction"){
        X <- matrix(rnorm(n=N*p,0,1), N, p)
        X[,p] <- X[,1]*X[,2]
    }

    if (regressors == "corr1"){
        mu <- rep(0,p)
        Sigma <- 1
        sigma <- matrix(rep(0.5),p,p)
        for (i in 1:p){
            sigma[i,i] <- Sigma
        }
        X <- mvrnorm(N, mu, sigma)
    }
    
    if (regressors == "corr2"){
        mu <- rep(0,p)
        Sigma <- 1
        sigma <- diag(Sigma, p, p)
        sigma[1,2] <- 0.9
        sigma[2,1] <- 0.9
        X <- mvrnorm(N, mu, sigma)
    }
    
    #### Treatment Probability
    
    W <- sample(c(0,1), N, replace = TRUE, prob = c(0.5, 0.5))
    
    if (regressors == "selectionbias"){
        X <- matrix(rnorm(n=N*p,0,1), N, p)
        W <- ifelse(X[,1] > 0, sample(c(0,1), N, replace = TRUE, prob = c(0.25, 0.75))
                ,sample(c(0,1), N, replace = TRUE, prob = c(0.75, 0.25)))
    }
    #### Beta Vector
    
    b_vec <- seq(0.5,1,length.out = p)
    
    if (regressors == "beta12"){
        X <- matrix(rnorm(n=N*p,0,1), N, p)
        b_vec[1:2] <- 0
    }
    
    #### true treatment effects 
    
    if (treatment == "treat1"){
        t <- 0.2
    }
    
    if (treatment == "treat2"){ 
        t <- 0.1 + 0.1*ifelse(X[,1] > 0, 1, 0) +  0.1*ifelse(X[,2] > 0, 1, 0)
    }
    
    if (treatment == "treat3"){ 
        t <- 0.1 + X[,1]*ifelse(X[,1] > 0, 1, 0) +  X[,2]*ifelse(X[,2] > 0, 1, 0)
    }
    
    if (treatment == "treat4"){ 
        t <- 0.1 + X[,1]*ifelse(X[,1] > 0, 1, 0) +  X[,2]*ifelse(X[,2] > 0, 1, 0) + X[,3]*ifelse(X[,3] > 0, 1, 0) + X[,4]*ifelse(X[,4] > 0, 1, 0)
    }
    
    if (treatment == "treat5"){
        t <- 0.2 + X[,1]*X[,2]
    }
    
    
    Y <- X %*% b_vec + t*W + rnorm(N,0,0.1)
    df <- data.frame(Y,W,t,X)
    
    return(df)
}




#-----------------------------------------------------------------------------
#### Causal Forest - Treatment Effects
#-----------------------------------------------------------------------------


cf_selectmedian <- function(df){
    Y <- df$Y
    W <- df$W
    X <- df[,-(1:3)]
    
    ### Untrainted CF
    cf_pilot <- causal_forest(X, Y, W)   #, num.trees = 5000)
    var_imp <- variable_importance(cf_pilot) 
    select_index <- which(var_imp >= median(var_imp))
    #cat("The most important Variables are", select_index, "\n")
    
    ### Trained CF
    cf <- causal_forest(X[, select_index], Y, W)
    tau_hat <- predict(cf)$predictions
    df["treat"] = tau_hat
    
    return(df)
}


#-----------------------------------------------------------------------------
#### Simulation Study (NEW) - Table (Comparison)
#-----------------------------------------------------------------------------


simulation <- function(p,numsim,regressors){
    
    MSE_treat <- matrix(rep(NaN),5,4)    
    
    #### Treatment 1 row
    x <- 0
    
    for (n in c(100,150,250,500)){
        
        x <- x + 1
        
        MSE_temp <- c() 
        
        for (i in 1:numsim){
    
            #### DGP
            data <- dgp(p,N=n,treatment="treat1",regressors)
        
            #### Get Treatment Values
            cfval <- cf_selectmedian(data)
        
            #### MSE for individual treatment effects
            MSE_temp[i] <- mean((cfval$t - cfval$treat)^2)
        }
        MSE_treat[1,x] <- mean(MSE_temp)
    }
    
    #### Treatment 2 row
    x <- 0
    
    for (n in c(100,150,250,500)){
        
        x <- x + 1
        
        MSE_temp <- c() 
        
        for (i in 1:numsim){
    
            #### DGP
            data <- dgp(p,N=n,treatment="treat2",regressors)
        
            #### Get Treatment Values
            cfval <- cf_selectmedian(data)
        
            #### MSE for individual treatment effects
            MSE_temp[i] <- mean((cfval$t - cfval$treat)^2)
        }
        MSE_treat[2,x] <- mean(MSE_temp)
    }
    
    #### Treatment 3 row
    x <- 0
    
    for (n in c(100,150,250,500)){
        
        x <- x + 1
        
        MSE_temp <- c() 
        
        for (i in 1:numsim){
    
            #### DGP
            data <- dgp(p,N=n,treatment="treat3",regressors)
        
            #### Get Treatment Values
            cfval <- cf_selectmedian(data)
        
            #### MSE for individual treatment effects
            MSE_temp[i] <- mean((cfval$t - cfval$treat)^2)
        }
        MSE_treat[3,x] <- mean(MSE_temp)
    }    
    
    #### Treatment 4 row
    x <- 0
    
    for (n in c(100,150,250,500)){
        
        x <- x + 1
        
        MSE_temp <- c() 
        
        for (i in 1:numsim){
    
            #### DGP
            data <- dgp(p,N=n,treatment="treat4",regressors)
        
            #### Get Treatment Values
            cfval <- cf_selectmedian(data)
        
            #### MSE for individual treatment effects
            MSE_temp[i] <- mean((cfval$t - cfval$treat)^2)
        }
        MSE_treat[4,x] <- mean(MSE_temp)
    }    
    
    #### Treatment 5 row
    x <- 0
    
    for (n in c(100,150,250,500)){
        
        x <- x + 1
        
        MSE_temp <- c() 
        
        for (i in 1:numsim){
    
            #### DGP
            data <- dgp(p,N=n,treatment="treat5",regressors)
        
            #### Get Treatment Values
            cfval <- cf_selectmedian(data)
        
            #### MSE for individual treatment effects
            MSE_temp[i] <- mean((cfval$t - cfval$treat)^2)
        }
        MSE_treat[5,x] <- mean(MSE_temp)
    }    
    
    MSE <- round(MSE_treat,3)
    
    output <- data.frame("Application" = c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4","Treatment 5"), 
                         "MSE" = MSE)
    colnames(output)[2:5] <- c("N = 100", "N = 150", "N = 250", "N = 500")
    rownames(output) <- NULL
    return(output)
}



#-----------------------------------------------------------------------------
#### Causal Forest - Graphical Treatment Heterogeneity
#-----------------------------------------------------------------------------

CF_graphics <- function(p,N,treat,regressors,numsim){
    
    #### containers for group_by means
    X1_0 <- c()
    X1_1 <- c()
    
    X2_0 <- c()
    X2_1 <- c()
    
    for (i in 1:numsim){
    
        #### DGP
        data <- dgp(p,N,treat,regressors)
        
        #### Causal Forest Values
        cfval <- cf_selectmedian(data)
        
        #### Save group_by treatment effect means
        X1_0[i] <- mean(cfval$treat[cfval$X1 <= 0])
        X1_1[i] <- mean(cfval$treat[cfval$X1 > 0])
        
        X2_0[i] <- mean(cfval$treat[cfval$X2 <= 0])
        X2_1[i] <- mean(cfval$treat[cfval$X1 > 0])
    }
    cf <- cbind(X1_0,X1_1,X2_0,X2_1)
    
    boxplot(cf[,1],cf[,2],cf[,3],cf[,4]
            #, xlab = "X1 < 0.5 | X1 > 0.5 | X2 < 0.5 | X2 > 0.5"
            ,names = c("X1 < 0","X1 > 0","X2 < 0","X2 > 0")
            , ylab = "Estimated CATE"
            , main = "Causal Forest - Monte Carlo Results")
}  


#-----------------------------------------------------------------------------
#### KNN - Graphical Treatment Heterogeneity
#-----------------------------------------------------------------------------

knn_simstudy <- function(k,p,N,numsim,dgp){
    
    #### containers for group_by means
    X1_0 <- c()
    X1_1 <- c()
    
    X2_0 <- c()
    X2_1 <- c()
    
    for (i in 1:numsim){
    
        #### DGP
        data <- dgp(p,N)
        
        #### KNN-values
        knnval <- get_knnval(k,data)
        
        #### Save group_by treatment effect means
        X1_0[i] <- mean(cfval$treat[cfval$X1 <= 0])
        X1_1[i] <- mean(cfval$treat[cfval$X1 > 0])
        
        X2_0[i] <- mean(cfval$treat[cfval$X2 <= 0])
        X2_1[i] <- mean(cfval$treat[cfval$X1 > 0])    
    }
    
    knn <- cbind(X1_0,X1_1,X2_0,X2_1)
    boxplot(knn[,1],knn[,2],knn[,3],knn[,4]
            #, xlab = "X1 < 0.5 | X1 > 0.5 | X2 < 0.5 | X2 > 0.5"
            ,names = c("X1 < 0","X1 > 0","X2 < 0","X2 > 0")
            , ylab = "Estimated CATE"
            , main = "KNN - Monte Carlo Results")
}