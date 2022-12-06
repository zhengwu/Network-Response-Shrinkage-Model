#utility functions

# ----------------------------------------------  
# Some functions to calculate MSE
f_mse <- function(beta){ 
  A_est  <- X_test%*%t(beta) 
  m      <- mean((A_test - A_est)^2)
  return(m)
}


f_mse2 <- function(beta){ 
  A_est  <- X_test%*%t(beta) 
  m      <- mean( ((A_test - A_est)^2)/(A_test^2) )
  return(m)
}

f_mse3 <- function(beta){ 
  A_est  <- X_test%*%t(beta) 
  m      <- mean( abs((A_test - A_est)) )
  return(m)
}

# Function to standardize A, get smaller A, get C and E.
A_standard <- function(X, U, V, top){
  
  C = matrix(0,V,V)
  for ( r in 2:V ){C[1:(r-1),r] = 1}
  C = which(C==1,TRUE)                    
  C = C[order(C[,1],C[,2]),]         
  
  A0  <- X%*%(U[,C[,1]]*U[,C[,2]])
  
  C0     <- C                       
  non0   <- colSums(A0)             
  C      <- C[which(non0 != 0),]    
  non01  <- non0[(non0 != 0)]      
  topsig <- tail(sort(non01), floor(sum(non01 != 0)*top))   
  E      <- C0[(non0 %in% topsig),] 
  E      <- rbind(E, E[,c(2,1)])
  
  A  <- A0[, which(non0 != 0)]      
  return(list(A, C, E))
}




# --------------------------- Function for Rank3 ---------------------------
A_standard3 <- function(X, U1, U2, U3, V, top){
  
  C = matrix(0,V,V)
  for ( r in 2:V ){C[1:(r-1),r] = 1}
  C = which(C==1,TRUE)                    
  C = C[order(C[,1],C[,2]),]         
  
  A0  <- X%*% ( U1[,C[,1]]*U1[,C[,2]]  +  U2[,C[,1]]*U2[,C[,2]]  + U3[,C[,1]]*U3[,C[,2]] )
  
  C0     <- C                       
  non0   <- colSums(A0)             
  C      <- C[which(non0 != 0),]    
  non01  <- non0[(non0 != 0)]      
  topsig <- tail(sort(non01), floor(sum(non01 != 0)*top))   
  E      <- C0[(non0 %in% topsig),]
  E      <- rbind(E, E[,c(2,1)])
  
  A  <- A0[, which(non0 != 0)]      
  return(list(A, C, E))
}

