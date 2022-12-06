# fit linear model using the first 12 covariates in X_0 (age, gender, 10 PCs)
# get the residual and use it as Y, and fit our model using all the SNPs

#############  data involved in this code;
# prepared_fa_networkdata.mat - obtained from "prepare_hcpbrainnetwork.m" (matlab script)
# id_covSNP.rdata - subjects with genetic information
# X0.rdata -  confundors, 12 covariates including age, gender and 10 PCs of SNP info 
# fake.SNP.rdata - faked SNP data; note that we can't share the real data due to data usage terms
# Bsize_new.rdata - encoding the grouping structure along genome


rm (list = ls())

library(R.matlab)
  FA  <- readMat("./HCP-YA/prepared_fa_networkdata.mat")

  A  <- cbind(FA$all.id, FA$all.network.vector)
  C  <- FA$no.zeros.pairs
  V  <- dim(FA$all.network)[1]
  id <- FA$all.id

  
# common id among covariates, SNP and network data
  load("./HCP-YA/id_covSNP.rdata")
  
  sum(id %in% id_covSNP)                          # 1010
  sum(id_covSNP %in% id)  
  
  sum(intersect(id_covSNP, id) == intersect(id, id_covSNP))     # 1010 
  id_common <- intersect(id_covSNP, id)
  #save(id_common, file = "id_common_FA.rdata") 


# save the data 
  A <- A[which(A[,1] %in% id_common),]
  net_fa  <- list(A, C, V)
  #save(net_fa, file = "net_fa.rdata")



###################################################
##################################################

#----------------------------- Use Original FA data 
# A and X covariates part 
  load("./HCP-YA/net_fa.rdata")
  load("./HCP-YA/X0.rdata")


# ----------------------- Preparation ---------------------------
# Network data
  A  <- net_fa[[1]]
  A  <- A[order(A[,1]), ]
  C  <- net_fa[[2]]
  W  <- nrow(C)
  V  <- net_fa[[3]]                             # dim(all_network)[1]


# covariates  
  X0  <- X0[order(X0$id), ]  
  X0  <- as.matrix(X0[complete.cases(X0),])
  print( paste0("Matched id for A and X: ",   nrow(A) == sum(A[,1] == X0[,1]) )  )   #check if the subject id are matched     
  X0  <- X0[, -1]
  X0  <- X0[, c(2,1,3:ncol(X0))]   # Switch gender and age columns




########################################################################
# --------------------- Step 1, linear regression ----------------------
# standardize X0
  X0  <- cbind(X0[,1], apply(X0[,2:ncol(X0)],  2, function(x) scale(x)))
  colnames(X0)[1] <- "gender"
  
  pred0 <- NULL
  for (pair in 2:ncol(A) ){
    d      <- as.data.frame(cbind(A[,pair], X0))
    colnames(d)[1] <- "y"
    m0     <- lm(y~ gender + age + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, d)
    pred0  <- cbind(pred0, predict(m0, d[,2:ncol(d)]))
    
    print(pair)
  }
  rownames(pred0) <- NULL


# Save the residuals   
  A_resid       <-cbind(A[,1], A[,-1] - pred0)
  net_fa_resid  <- list(A_resid, C, V)
  #save(net_fa_resid, file = "net_fa_resid.rdata")



  
  
####################################################################  
# --------------------- Step 2, fit our model ----------------------
  
#-------------------------------------
  times <-  10            # random seed
  options(digits = 12) 
#-------------------------------------  
  
# A (residual after fitting linear models)
  load("./HCP-YA/net_fa_resid.rdata")
  
# SNP part  
  load("./HCP-YA/fake.SNP.rdata")
  
  
  load("./HCP-YA/Bsize_new.rdata") 
# Last one in Bsize_new is 1035, break the last Bsize into 50 SNP each.
  lastone <- Bsize_new[length(Bsize_new)]
  if((lastone %% 50) != 0 ){
    newsize <- c(rep(50, lastone/50), lastone %% 50)
  }else{
    newsize <- rep(50, lastone/50)
  } 
  Bsize_new <- c(Bsize_new[-length(Bsize_new)], newsize)
  
  
# load prior
 # library(R.matlab)
 #  E  <- readMat("network_prior.mat")
  E  <- FA$E.no.zeros.pairs
  
  
# ----------------------- Preparation ---------------------------
# Network data
  A  <- net_fa_resid[[1]]
  A  <- A[order(A[,1]), ]
  C  <- net_fa_resid[[2]]
  W  <- nrow(C)
  V  <- net_fa_resid[[3]]                         # dim(all_network)[1]
  

# SNPs     
  X1  <- cbind(SNP$id, SNP[, -ncol(SNP)])
  X1  <- as.matrix(X1[complete.cases(X1),])
  X1  <- X1[order(X1[, 1]), ]
  print( paste0("Matched id for A and SNP: ", nrow(A) == sum(A[,1] == X1[,1]) )  ) 
  X1  <- X1[, -1]
  
  
# Choose the training, test and validation set
  set.seed(times)
  print(paste0("seed: ", times))
  id      <- sample(A[,1], replace = F)
  
  ytrain  <- A[A[,1] %in% id[1:506], -1]
  ytest   <- A[A[,1] %in% id[507:758], -1]
  yvalid  <- A[A[,1] %in% id[759:length(id)], -1]   
  subid   <- A[A[,1] %in% id[759:length(id)], 1]
  
  xtrain  <- X1[A[,1] %in% id[1:506], ]
  xtest   <- X1[A[,1] %in% id[507:758], ]
  xvalid  <- X1[A[,1] %in% id[759:length(id)], ]
  
# Standardize (Use original Y, standardized X) 
  xtrain <- apply(xtrain, 2, function(x) scale(x)) 
  xtest  <- apply(xtest,  2, function(x) scale(x))
  xvalid <- apply(xvalid, 2, function(x) scale(x)) 
  
  print(c(dim(ytrain), dim(ytest), dim(yvalid),
          dim(xtrain), dim(xtest), dim(xvalid) ))  
  
# ----------------------
  N    <- nrow(xtrain)                          # sample size
  P    <- ncol(xtrain)                          # parameter size
  Q    <- length(Bsize_new)
  size <- Bsize_new
  PW   <- rep(1:Q, size)  
  
# ---------------------------------------------------------------------------    
  

  
# Fit the model for different combination of parameters
# ------------------------------- Fit Model -------------------------------------
  source("BConNet_theta.r") 
  h  = 1
  f  = 2
  mu = 1
    
    
  lam_H = matrix(rep(h, ncol(X1)),  nrow = P,  ncol = V)
  lam_f = rep(f,  length(Bsize_new)) 
  nu    = 0.1
  print(paste0("para: ", c(h, f, mu)))
  
  res   = BConNet(ytrain, C, xtrain, mu,nu,lam_f,lam_H,PW,E,a_sigma=1,b_sigma=1,a_omega=2,b_omega=1, inititer=50)
  out   <- res$U
  theta <- res$theta  
  
# Save U and theta for test (valid) set 
  est     <- out[,C[,1]]*out[,C[,2]]
  theta1  <- matrix(rep(theta, each = nrow(xtest)), nrow = nrow(xtest))
  A_est   <- theta1 + xtest %*% est 
  mse     <- mean((ytest - A_est)^2)
  bwise   <- sum(rowSums(abs(out)) != 0)
  corr    <- cor(c(ytest), c(A_est))
  
# Calculate the correlation between rows in ytest and estimation
  corr_mean <- NULL
  for(row in 1:nrow(A_est)){  corr_mean <- c(corr_mean, cor(ytest[row,], A_est[row,])) }
  corr_mean <- mean(corr_mean)   
  
  h1      <- h
  f1      <- f
  U       <- c(list(out), list(theta))
  
  
  save(U, file = paste0("Utrain_", 
                        h1, "_", f1, "_", mu, "_", times, ".rdata"))
  
# Save ytest and estimation by xtest  
  test <- c(list(ytest), list(A_est))
  save(test, file = paste0("test_", 
                           h1, "_", f1, "_", mu, "_", times, ".rdata"))
  
  
###########################################################################
#  Validaion
  
# calculate MSE for valid (test) set  
  est     <- out[,C[,1]]*out[,C[,2]]
  theta1  <- matrix(rep(theta, each = nrow(xvalid)), nrow = nrow(xvalid))
  A_est   <- theta1 + xvalid %*% est 
  mse1    <- mean((yvalid - A_est)^2)
  
  corr1        <- cor(c(yvalid), c(A_est))
  corr1_mean   <- NULL
  for(row in 1:nrow(A_est)){  corr1_mean <- c(corr1_mean, cor(yvalid[row,], A_est[row,])) }
  corr1_mean <- mean(corr1_mean)
  
# save yvalid and A_est
  y <- c(list(subid), list(yvalid), list(A_est))
  save(y, file = paste0("y_", 
                        h1, "_", f1, "_", mu, "_", times, ".rdata"))
  
  
# ------------------------------------------------
#calculate Correlation between Y and Y_estimate
  MSE   <- c(h1, f1, mu, 
             mse, bwise, mse1, 
             corr, corr_mean, corr1, corr1_mean, times)
  save(MSE, file = paste0("MSE_", 
                          h1, "_", f1, "_", mu, "_", times, ".rdata"))
  print(MSE)
  
  
  
  
  