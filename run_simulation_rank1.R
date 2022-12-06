

#############################################################
#                       Run Models                          #         
#############################################################  
rm(list = ls())
source('utility_functions.R')

#Run on simulated data (Rank 1)

######################################################################################  
#load U (which is scriptB in formula (1) in the manuscript)

# ------------------------------ Rank 1 ----------------------------------------------
dset   <-   1  # index for the dth dataset  (any number between 1-150)
V_size <-   320  # 320: U3 V=20       370: U3 V=70        420: U4 V=20        470: U4 V=70
# ------------------------------------------------------------------------------------    

if(V_size == 320){
  load("./rank1/U3_V20_random.RData")
  V = 20
  U = t(U)
  
}else if(V_size == 370){
  load("./rank1/U3_V70_random.RData")
  V = 70
  U = t(U)
  
}else if(V_size == 420){
  load("./rank1/U4_V20_random.RData")
  V = 20
  U = t(U)
  
}else if(V_size == 470){
  load("./rank1/U4_V70_random.RData")
  V = 70
  U = t(U)
  
}else{
  print("no V size")}   


load("SNP_all.rdata")
X       <- SNP_all[[dset]]
X_test  <- SNP_all[[dset+150]]


# -----------------------------------------------------------------------------
# ----------------------- Variable that will change ---------------------------
set.seed(dset)
sigma2 <- 0.1                                # 0.1 or 0.5 （for noise level)
top    <- 0.7                             # percentage of columns that we want to use as prior

N  = 100
P  = 2000
Q  = P/100
PW = rep(1:Q,each=100)

out   <- A_standard(X, U, V, top)
C     <- out[[2]]
W     <- nrow(C)


Error <- rnorm(N*W, 0, sqrt(sigma2));   print(sum(Error))
A     <- out[[1]] + Error
# ------------------------------------------------------------------------------



######################################### 0. Run Our Model #################################################
source("BConNet.r")  

lam_H  =   1
lam_f  =   300
mu     =  5
nu     = 0.1


# without prior
E = NULL

t_est_0 <- system.time({  
  res_0  = BConNet(A,C,X,mu,nu,lam_f,lam_H,PW,E,a_sigma=1,b_sigma=1,a_omega=2,b_omega=1)
})
out0  <- res_0$U
est0  <- out0[,C[,1]]*out0[,C[,2]]
est0  <- t(est0)


# with prior
E = out[[3]] 

t_est_1 <- system.time({ 
  res_1  = BConNet(A,C,X,mu,nu,lam_f,lam_H,PW,E,a_sigma=1,b_sigma=1,a_omega=2,b_omega=1)
})
out1  <- res_1$U
est1  <- out1[,C[,1]]*out1[,C[,2]]
est1  <- t(est1)


###############################################  1. SGL  ####################################################
library(SGL)

t_1  <- system.time({
  
  index  <- PW                          
  B1     <- NULL
  for (pair in 1:ncol(A)){
    y    <- A[,pair]
    data <- list(x = X, y = y)
    fit1 <- SGL(data, index, type = "linear", lambdas = 0.0045)
    B1   <- rbind(B1, fit1$beta)
    print(pair)
  }
  
})  


#############################################   2. EMSHS   ###################################################  
library(EMSHS)

t_2 <- system.time({      
  
  B2     <- NULL
  for (pair in 1:ncol(A)){
    y    <- A[,pair]
    fit2 <- EMSHS(y, X, 5, nu, E = E)         # use mu = 5
    B2   <- rbind(B2, t(fit2$beta))
    print(pair)
  }
  
})



############################################# 4. Blasso #####################################################
require(monomvn)

t_4 <-  system.time({ 
  B4     <- NULL 
  for (pair in 1:ncol(A)){
    y    <- A[,pair]
    fit4 <- blasso(X = X, y = y, T = 200)  
    aic  <- 2*apply(fit4$beta, 1, function(x) sum(x != 0)) -  2*fit4$llik     
    B4   <- rbind(B4, fit4$beta[which.min(aic),]) 
    print(paste0("min AIC id: ", which.min(aic)))
    print(pair)
  }
  
}) 


# -----------------------------------------------------------------------------------------------------------  
# Test

out    <- A_standard(X_test, U, V, top)
C      <- out[[2]]
E      <- out[[3]]
W      <- nrow(C)

dset1  <- dset + 200       
set.seed(dset1)
Error1 <- rnorm(N*W, 0, sqrt(sigma2));    print(sum(Error1))
A_test <- out[[1]] + Error1 


origin <- U[,C[,1]]*U[,C[,2]]
origin <- t(origin)
# save beta to project file 
print(c( dim(origin), dim(est1), dim(B1), dim(B2), dim(B4) ) )
b <- list(origin, est0, est1, B1, B2, B4)

# Time
t  <- rbind(t_est_0, t_est_1, t_1, t_2, t_4)[, c(1,3)]

if(V_size %in% c(320, 370, 420, 470)){
  save(b, file = paste0("first8_beta", V_size, "_", sigma2, "_", dset, ".rdata"))
  save(t, file = paste0("first8_t", V_size, "_", sigma2, "_", dset, ".rdata"))
}else{
  save(b, file = paste0("last8_beta", V_size, "_", sigma2, "_", dset, ".rdata"))
  save(t, file = paste0("last8_t", V_size, "_", sigma2, "_", dset, ".rdata"))
}



# Calculate MSE for beta
d1 <- cbind(c(origin), c(est0), c(est1), c(B1), c(B2), c(B4))
mb <- apply(d1, 2, function(x) mean((d1[,1] - x)^2))[-1]      


# Calculate prediction MSE
m     <- c(f_mse(est0), f_mse(est1), f_mse(B1), f_mse(B2), f_mse(B4))   
m_new <- c(f_mse2(est0), f_mse2(est1), f_mse2(B1), f_mse2(B2), f_mse2(B4))      
m_abs <- c(f_mse3(est0), f_mse3(est1), f_mse3(B1), f_mse3(B2), f_mse3(B4))          

# Calculate AUC
library(pROC)
d1[,1] <- ifelse(d1[,1] != 0, 1, 0)
AUC    <- apply(d1, 2, function(x)  auc(roc(d1[,1], abs(x))) )[-1]           

# Data set number  & V_size  & sigma2
dindex  <- rep(dset, 5)  
setting <- rep(V_size, 5)
sig2    <- rep(sigma2, 5)

# combine
model       <- c("OurModel(noPrior)", "OurModel(Prior)", "SGL", "EMSHS", "Blasso")
a           <- cbind(model, t, mb, m, m_new, m_abs, AUC, dindex, setting, sig2)
colnames(a) <- c("model", "time(user)", "elapsed", "MSE(beta)", "MSE(prediction)", "MSE(new)", "MSE(abs)", "AUC", "dataset", "setting", "sigma2")

print(a)
#save(a, file = paste0("a", V_size, "_", sigma2, "_", dset, ".rdata")) 
