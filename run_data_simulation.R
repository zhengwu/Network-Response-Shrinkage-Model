####################################################
#                 generate data                    #
####################################################

library(glmnet)
library(ncvreg)
library(R.matlab)
library(SGL)


load("befgen4.RData")


# combining every 3 SNP size
  l <- seq(1, length(LD_idx)+1, 3)
  LD_idx_new <- NULL
  Bsize_new  <- NULL
  for (i in l) {
    if(i == length(LD_idx)){
      LD_idx_new <- c(LD_idx_new, LD_idx[i])
    }
    li <- list(unlist(c(LD_idx[i], LD_idx[i+1], LD_idx[i+2])))     
    Bsize_new <- c(Bsize_new, length(li[[1]]))                                          
    LD_idx_new <- c(LD_idx_new, li)
  }

# combining every SNP size until the size exceeds 100
  LD_idx_new1  <- NULL
  Bsize_new1   <- NULL

  i = 1                 # indicate the ith list in LD_idx_new
  j = 0                 # indicate the jth list in LD_idx_new1
  repeat{
    j = j+1       
    Bsize_new1 <- c(Bsize_new1, length(LD_idx_new[[i]]))
    list_i     <- LD_idx_new[i]
    
    
    if(Bsize_new1[j] < 100){ 
      while(Bsize_new1[j] < 100){
        li     <- list(unlist(c(list_i, LD_idx_new[i+1])))
        list_i <- li
        Bsize_new1[j] <- Bsize_new1[j] + length(LD_idx_new[[i+1]])
        i      <- i + 1
        print(i)
        if(i == length(LD_idx_new)){break}
      }
    }
    
    i <- i + 1
    LD_idx_new1 <- c(LD_idx_new1, list_i)
    if(i >= length(LD_idx_new)){break}
  }



############################### Generate data ###########################################
#dset <- as.numeric(args[1]) 
dset = 1 
set.seed(dset)

# -------------------------------------------------------------------
  N=ncol(GT2hap)  #sample size used to generate LD block
  N1 = 100        #sample size:   100training + 100testing
  NB1 = 50        #Number of groups 
  NB2 = 100       #Number of snps in each group

# ---------------------------------------------------------------------

  NS1=NB2*NB1
  
  idx0=which(Bsize_new1>= 100)
  idx=sort(sample(idx0,NB1,replace=F))    # sample NB1(=50) LD blocks from all blocks that has block size >= 100
  
  
  Bsizeo=Bsize_new1[idx]
  LD_idx_o=LD_idx_new1[idx]
  LD_idx_o_e=LD_idx_o
  
  for(i in 1:NB1){                       
    idx1=sort(sample.int(Bsizeo[i],NB2))
    
    LD_idx_o_e[[i]]=LD_idx_o[[i]][idx1]
    print(length(LD_idx_o_e[[i]]))
    
    print(i)
  }
  
  
  col_sel_e=unlist(LD_idx_o_e)
  
  NS=length(col_sel_e)
  
  SNP=matrix(rep(0,N1*NS),ncol=NS)
  

# --------------------------------------------------------------- 
  k=0
  for(i in 1:NB1){
    idx2 = sample.int(N,N1*2,replace=T)       
    SNP[,(k+1):(k+100)]=t(GT2hap[LD_idx_o_e[[i]],idx2[1:N1]]+        
                            GT2hap[LD_idx_o_e[[i]],idx2[1:N1+N1]])     
    k=k+100
  }

  
# --------------------------------------------------------------
  X <- SNP
  
  p  = ncol(X)
  Id = matrix(0,p,p)
  rc = NULL
  for ( i in 1:(p-1) ){
    for ( j in (i+1):p ){
      if ( sum(X[,i]!=X[,j]) == 0 ){
        Id[i,j] = 1
        rc  = rbind(rc, c(i,j))
      }} 
  }
  sum(Id)
  
  # find the column number for the repeated column
  repcol <- unique(c(rc[,1],rc[,2]))
  SNP    <- X[, -repcol]
  SNP    <- SNP[,1:2000]
  
  # change the group index for the gene so that we have 20 groups in total
  group_id  <- rep(1:NB1,each=NB2)
  
  group_id1 <- group_id[-repcol][1:2000]
  need      <- 100 - table(group_id1[which(group_id1 <= 20)])    
  PW        <- c(group_id1[which(group_id1 <= 20)], rep(1:20, need))
  
  
  #save(SNP, file = paste0("SNP_", dset, ".rdata"))
  #save(PW,  file = paste0("PW_", dset, ".rdata"))
  
  
#############################################################
#                       Generate U                          #         
#############################################################
# Note:  V = 20 or 70 should use different seed, otherwise the number of 1 in each SNP will be the same
# seed1000 for V=20, seed100 for V=70   
  
s = 1000 
#s = 100  
set.seed(s)
V = 20
#V = 70                          
#------------------------------
  P  = 2000
  Q  = P/100
  PW = rep(1:Q,each=100)
  
  C  = matrix(0,V,V)
  for ( r in 2:V ) {C[1:(r-1),r] = 1}
  
  
  C  = which(C==1,TRUE)                      
  C  = C[order(C[,1],C[,2]),]         
  W  = nrow(C)
  
###################################### Random Column in each block #########################################
  loc <- function(non0_id){
    a <- if(non0_id == 20) 1:20   else  ((idex - floor(non0_id/2)):(idex + floor(non0_id/2)))[1:non0_id]
    return(a)}
  
# ----------------------------------------- Signal Type I -----------------------------------------------
  U      = matrix(0,V,P)
  idex   = ceiling(nrow(U)/2) 
  sig    = 100                 
  
  if(V == 20){
    non0 = round(rnorm(sig, 15, 1.9))
  }else{non0 = round(rnorm(sig, 35, 1.9))}         # for V = 20, sample from N(15, 1.9), V = 70 sample from N(35, 1.9)
  snp_id = replicate(2, sample(1:100, 50, FALSE))
  
# ---------------------------------------- 
# Signal all = 1  
  for(snp in 1:(sig/2)){       
    col_id <- snp_id[snp,]
    
    row  = loc(non0[snp])
    row1 = loc(non0[snp + 50])
    
    U[row,  col_id[1]]      = rep(1, non0[snp])
    U[row1, col_id[2]+100]  = rep(1, non0[snp + 50])
  }
  # check
  a  <- colSums(U)
  a1 <- t(matrix(a, ncol = 20))
  apply(a1, 1, function(x) sum(x != 0))
  
save(U, file = paste0("U1_V", V,"_random.RData"))  
  
  
# ---------------------------------------- 
# Signal from N(0, 0.5^2)
  #U = matrix(0,V,P)
  if(V == 20){
    load("U1_V20_random.RData")
  }else{
    load("U1_V70_random.RData")
  }
  U1 <- U
  U  <- matrix(0,V,P)
  
  non0   <- apply(U1, 2, function(x) sum(x !=0))[1:200]
  non0   <- non0[which(non0 !=0)]
  snp_id <- matrix(which(apply(U1, 2, function(x) sum(x !=0))[1:200] != 0), ncol = 2)
# ----------------------
  sigma2 = 0.5
# ----------------------
  
  for(snp in 1:(sig/2)){       
    col_id <- snp_id[snp,]
    
    row  = loc(non0[snp])
    row1 = loc(non0[snp + 50])
    
    U[row,  col_id[1]]  = rnorm(non0[snp],      mean = 0, sd = sqrt(sigma2))
    U[row1, col_id[2]]  = rnorm(non0[snp + 50], mean = 0, sd = sqrt(sigma2))
    
  }
  # check
  a <- colSums(U)
  a1 <- t(matrix(a, ncol = 20))
  apply(a1, 1, function(x) sum(x != 0))
  sd(c(U[U!=0]))^2  
  
save(U, file = paste0("U3_V", V,"_random.RData"))
  
  
# ---------------------------------- Signal Type III (Random SNP position) ---------------------------------------
  U      = matrix(0,V,P)
  idex   = ceiling(nrow(U)/2)
  sig    = 100
  
  if(V == 20){
    non0 = round(rnorm(sig, 15, 1.9))
  }else{non0 = round(rnorm(sig, 35, 1.9))}    
  
  snp_id = replicate(5, sample(1:100, 20, FALSE))     
  
  
# ---------------------------------- Signal all = 1      
  for(snp in 1:(sig/5)){
    col_id <- snp_id[snp,]
    
    row1 = loc(non0[snp]);    row2  = loc(non0[snp+20])
    row3 = loc(non0[snp+40]); row4  = loc(non0[snp+60])
    row5 = loc(non0[snp+80])
    
    U[row1, col_id[1]]     = rep(1, non0[snp]);     U[row2, col_id[2]+100]  = rep(1, non0[snp+20])
    U[row3, col_id[3]+200] = rep(1, non0[snp+40]);  U[row4, col_id[4]+300]  = rep(1, non0[snp+60])
    U[row5, col_id[5]+400] = rep(1, non0[snp+80])
    
  }
  
# check 
  a  <- colSums(U)
  a1 <- t(matrix(a, ncol = 20))
  apply(a1, 1, function(x) sum(x != 0))
  
save(U, file = paste0("U2_V", V,"_random.RData"))
  
  
  
# -------------------------------- Signal from N(0, 0.5)
  if(V == 20){
    load("U2_V20_random.RData")
  }else{
    load("U2_V70_random.RData")
  }
  U1 <- U
  U  <- matrix(0,V,P)
  
  non0   <- apply(U1, 2, function(x) sum(x !=0))[1:500]
  non0   <- non0[which(non0 !=0)]
  snp_id <- matrix(which(apply(U1, 2, function(x) sum(x !=0))[1:500] != 0), ncol = 5)
# -----------------------  
  sigma2    <- 0.5
# -----------------------     
  
  for(snp in 1:(sig/5)){
    col_id <- snp_id[snp,]
    
    row1 = loc(non0[snp]);    row2  = loc(non0[snp+20])
    row3 = loc(non0[snp+40]); row4  = loc(non0[snp+60])
    row5 = loc(non0[snp+80])
    
    U[row1, col_id[1]] = rnorm(non0[snp],    mean = 0, sd = sqrt(sigma2));       U[row2, col_id[2]]  = rnorm(non0[snp+20], mean = 0, sd = sqrt(sigma2))
    U[row3, col_id[3]] = rnorm(non0[snp+40], mean = 0, sd = sqrt(sigma2));       U[row4, col_id[4]]  = rnorm(non0[snp+60], mean = 0, sd = sqrt(sigma2))
    U[row5, col_id[5]] = rnorm(non0[snp+80], mean = 0, sd = sqrt(sigma2))
    
  }
  
# check 
  a  <- colSums(U)
  a1 <- t(matrix(a, ncol = 20))
  apply(a1, 1, function(x) sum(x != 0))
# check if the nonzero place are the same  
  # sum(apply(U, 2, function(x) sum(x !=0))[1:500] == apply(U1, 2, function(x) sum(x !=0))[1:500])   
  # sd(c(U[U!=0]))^2
  # 
  # save(U, file = paste0("U4_V", V,"_random.RData"))  
  # 
  # 
  # load("U4_V20_random.RData") 
  #   U1 <- U
  # load("U4_V20_rank3_2.RData")
  #   U2 <- U
  # load("U4_V20_rank3_3.RData") 
  #   U3 <- U
  #   
  # U_rank3 <- list(U1, U2, U3)
  # 
  # save(U_rank3, file = "U4_V20_rank3.rdata")
  # 
  

  