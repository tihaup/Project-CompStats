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
#### KNN - Treatment Effects
#-----------------------------------------------------------------------------

get_knnval <- function(k,df){
    
    #### Importing Data
    
    #### Brauchen wir df_train und df_test????
    df_train <- df
    rownames(df_train) <- 1:nrow(df_train)
    
    #### Splitting into Treated and Untreated Training Dataset
    df_t <- df_train[which(df_train$W == 1),]
    df_c <- df_train[which(df_train$W == 0),]

    drop_cols <- c("Y", "W", "t")
    
    #### calculate indices of k treated/untreated nearest neighbours
    treated_nn <- get.knnx(df_t[, !(names(df_t) %in% drop_cols)],                
                          query = df_train[, !(names(df_train) %in% drop_cols)], 
                          k = k)
    treated_nn <- data.frame(treated_nn$nn.index)

    untreated_nn <- get.knnx(df_c[, !(names(df_c) %in% drop_cols)],              
                          query = df_train[, !(names(df_train) %in% drop_cols)], 
                          k = k)
    untreated_nn <- data.frame(untreated_nn$nn.index)

    #### replacing index values
    treated_nn <- sapply(treated_nn, FUN = function(x){
                        df_t$Y[x]
                })
    
    untreated_nn <- sapply(untreated_nn, FUN = function(x){
                        df_c$Y[x]
                })


    #### Save "Potential Outcomes W/o Treatment in dataset": Y1 and Y0
    for (i in 1:nrow(df_train)){
        df_train$Y1[i] <- mean(treated_nn[i,])
        df_train$Y0[i] <- mean(untreated_nn[i,])
        df_train$treat <- df_train$Y1 - df_train$Y0
    }
    
    return(df_train)

}


#-----------------------------------------------------------------------------
#### Simulation Study - Table (Comparison)
#-----------------------------------------------------------------------------

simstudy <- function(p,numsim,treatment,regressors){
    
    MSE_cf <- c()
    MSE_10nn <- c()
    MSE_20nn <- c()    
    
    x <- 0
    
    for (n in c(100,150,250,500)){
        
        x <- x + 1
        
        MSE1_cf <- c()
        MSE1_10nn <- c()
        MSE1_20nn <- c()    
        
        for (i in 1:numsim){
    
            #### DGP
            data <- dgp(p,N=n,treatment,regressors)
        
            #### Get Treatment Values
            cfval <- cf_selectmedian(data)
            knnval10 <- get_knnval(10,data)
            knnval20 <- get_knnval(20,data)
        
            #### MSE for individual treatment effects
            MSE1_cf[i] <- mean((cfval$t - cfval$treat)^2)
            MSE1_10nn[i] <- mean((knnval10$t - knnval10$treat)^2)
            MSE1_20nn[i] <- mean((knnval20$t - knnval20$treat)^2)
        }
        MSE_cf[x] <- mean(MSE1_cf)
        MSE_10nn[x] <- mean(MSE1_10nn)
        MSE_20nn[x] <- mean(MSE1_20nn)
    }
    MSE <- cbind(MSE_cf, MSE_10nn, MSE_20nn)
    
    output <- data.frame(output <- data.frame("Method" = c("Causal Forest", "KNN-10", "KNN-20"), 
                                              "100" = c(round(MSE[1,],3)),
                                              "X150" = c(round(MSE[2,],3)),
                                              "X250" = c(round(MSE[3,],3)),
                                              "X500" = c(round(MSE[4,],3))
                                             ))
    colnames(output)[2:5] <- c("N = 100", "N = 150", "N = 250", "N = 500")
    rownames(output) <- NULL
    return(output)
}



#-----------------------------------------------------------------------------
#### Causal Forest - Graphical Treatment Heterogeneity
#-----------------------------------------------------------------------------

CF_simstudy <- function(p,N,numsim,dgp){
    
    #### containers for group_by means
    X1_0 <- c()
    X1_1 <- c()
    
    X2_0 <- c()
    X2_1 <- c()
    
    for (i in 1:numsim){
    
        #### DGP
        data <- dgp(p,N)
        
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