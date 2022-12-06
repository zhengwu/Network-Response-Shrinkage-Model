# Function BConNet
# written by Changgee Chang
# ver. 20190918
#
# A:       N by W brain network connectivity
# C:       W by 2 matrix connected pairs of ROI
# X:       N by P predictors
# mu:      ROI-wise shrinkage paramenter
# nu:      ROI-wise shrinkage adaptivity parameter
# lam_f:   SNP-wise penalty
# lam_H:   penalty for individual coefficients
# PW:      P by 1 vector indicating the group that each gene belongs to (1,2,...,Q)
# E:       e by 2 matrix with edges, indicating a graph between ROIs.
# a_sigma, b_sigma, a_omega, b_omega: prior parameters
# eps:     algorithm stops if relative improvement goes below eps


BConNet <- function(A,C,X,mu,nu,lam_f,lam_H,PW,E=NULL,a_sigma=1,b_sigma=1,a_omega=2,b_omega=1,inititer=25,gradual=5,eps=1e-4)
{
  N = nrow(A)
  W = ncol(A)                              #W: number of pairs of brain region
  
  if ( nrow(C) != W )
    stop("Invalid dimension of C")
  
  V = max(C)
  C2 = rbind(C,C[,2:1])
  C2 = C2[order(C2[,1],C2[,2]),]
  
  # diff() default:   calculate the difference between two numbers next to each other,   X_{i+1} - X_{i}
    Cidx = c(0,which(diff(C2[,1])!=0),2*W)
    Cmap = matrix(0,V,V)
    Cmap[C] = 1:W                            # assign value from 1 to W to the position (upper part of Cmap) indicated in C  
    Cmap = Cmap + t(Cmap)                    # make Cmap symmetric
    
    Cmat = matrix(0,2*W,V)
    Cmat[cbind(1:(2*W),C2[,1])] = 1
  
  
# ---------------------------------------------------------------------------------------------------
  if ( nrow(X) != N )
    stop("Dimensions of A and X do not match")
  P = ncol(X)                                 # number of SNP
  
  if ( length(PW) != P )
    stop("Dimensions of X and PW do not match")
  Q = max(PW)                                 # number of groups of SNP
  PWidx = c(0,which(diff(PW)!=0),P)

  # Graph edges
  if ( is.null(E) )
    e = 0
  else
  {
    e = nrow(E)
    E = E[order(E[,1],E[,2]),]
  }

  
 # ---------------------------------------------------------------------------------------------------
 # Pre-computations
   tX   = t(X)
   Asum = apply(A,2,sum)                     # column sum of NxW A matrix   (col sum for each pair of brain region)
   Abar = Asum/N
   Xsum = apply(X,2,sum)                     # column sum for each SNP 
   Xbar = Xsum/N 
   
   AA   = apply(A^2,2,sum)                   # column sum of square for each pair of brain region
   XA   = tX%*%A                             # (P x N) x (N x W)   =  P x W       SNP and Pair of brain region
  
  
 # Initialize f and g
   f     = rep(1,Q)                          # create f for each pair of brain region (Q pairs in total)
   sqf   = sqrt(f)
   alpha = rep(0,V)
   lam   = exp(alpha)                        # ?????????
   
   g     = rep(1,V)                          # g_{r} is for every brain region, there are V brain regions 
  
 # Initialize H
   H = matrix(0.01,P,V)                      # h_{rp} for each combination of brain region and SNP

  
   gH    = H*rep(g,each=P)                   # repeat each element in g (one by one) P times then multiply 
                                             # them with each element in H (col by col). 
                                             # ex: g = (g1, g2, ...., gv),   
                                             #     rep(g, each = P) = (g1_1, g1_2, ..., g1_p, ...., gv_1, gv_2,.., gv_p)
                                             # then H*rep(g, each = P) means each element in col1 is multiplied by g1
                                             # each element in col2 is multiplied by g2,... each element in colV is multiplied by gv
   U     = gH*sqf[PW]                        # Matrix form for u_{rp}
   
   #-----------------------------------------------------------------
   gHrs  = gH[,C[,1]]*gH[,C[,2]]             #????????????????
   Urs   = gHrs*f[PW]
   
   XUrs  = X%*%Urs
   XXUrs = tX%*%XUrs

   gr_f = rep(0,Q)
   gr_g = rep(0,V)
   gr_H = matrix(0,P,V)
  
   niter  = 0
   step_H = 0.001
   step_g = 0.001
   step_f = 0.001
  ###############################################################################
   repeat
   {
    niter = niter + 1
    pU = U
    
 # ----------------------------- E-Steps ----------------------------------------
 # E-Step: Omega
    Eomega = 2*nu*a_omega / (2*nu*b_omega + (alpha[E[,1]]-alpha[E[,2]])^2)   #page4 
    EOmega = matrix(0,V,V)
    EOmega[E] = -Eomega
    diag(EOmega) = 1 - apply(EOmega,1,sum)
    
    
    
 # ----------------------------- M-Step ------------------------------------------
 # M-Step: alpha
    if ( niter <= inititer )
      grmu = 0
    else
      grmu = mu*min(1,(niter-inititer)/gradual)
    
    H_alpha  = EOmega/nu + diag(lam*abs(g),V)
    tmp      = drop(EOmega%*%(alpha-grmu))/nu
    gr_alpha = tmp - rep(1,V) + lam*abs(g)
    L        = sum(tmp*(alpha-grmu))/2 - sum(alpha) + sum(lam*abs(g)) 
    
    cH     = chol(H_alpha)
    dir    = backsolve(cH,forwardsolve(t(cH),gr_alpha))
    maxdir = max(abs(dir))
    m      = sum(gr_alpha*dir)
    
    ss = 1
    repeat
    {
      nalpha = alpha - ss*dir
      nlam   = exp(nalpha)
      tmp    = drop(EOmega%*%(nalpha-grmu))/nu
      nL     = sum(tmp*(nalpha-grmu))/2 - sum(nalpha) + sum(nlam*abs(g))
      
      if ( ss*maxdir < eps | L-nL > ss*m*0.05 )
      {
        alpha = nalpha
        lam   = nlam
        break
      }
      ss = ss/2
    }      
    
    
  # --------------------------------------------------------------------------------------------  
  # Update theta and sig2
    UrsXbar = Urs*Xbar
    UrsXA   = Urs*XA
    theta   = Abar - apply(UrsXbar,2,sum)
    R2      = AA - 2*apply(UrsXA,2,sum) + apply(XUrs^2,2,sum) - N*theta^2
    sig2    = (sum(R2)+b_sigma)/(N*W+a_sigma)
    Xtheta  = outer(Xsum,theta)
    
  # --------------------------------------------------------------------------------------------  
  # Update H
    
  # Gradients
    gr_H = ((XXUrs-XA+Xtheta)[,Cmap[C2]]*U[,C2[,2]]) %*% Cmat
    gr_H = gr_H*outer(sqf[PW],g)

  # Proximal gradient descent with backtracking
    L = sum((A-rep(theta,each=N)-XUrs)^2)
    repeat
    {
      if ( step_H < 1e-7 )
        break
      
      tmpH = H - step_H * gr_H
      if ( niter <= inititer )
        nH = sign(tmpH)*pmax(abs(tmpH)-step_H*lam_H/100000,0)
      else
        nH = sign(tmpH)*pmax(abs(tmpH)-step_H*lam_H*min(1,(niter-inititer)/gradual),0)

      G_H = H - nH
      
      ngH = nH*rep(g,each=P)
      nU  = ngH*sqf[PW]
      
      ngHrs = ngH[,C[,1]]*ngH[,C[,2]]
      nUrs  = ngHrs*f[PW]
      nXUrs = X%*%nUrs
      nL    = sum((A-rep(theta,each=N)-nXUrs)^2)
      
      term2 = sum(gr_H*G_H)
      term3 = sum(G_H*G_H)
      if ( nL < L - term2 + term3/step_H/2 )
      {
        H  = nH
        gH = ngH
        U  = nU
        gHrs  = ngHrs
        Urs   = nUrs
        XUrs  = nXUrs
        XXUrs = tX%*%XUrs
        break
      }
      step_H = step_H/2
    }
    step_H = step_H*1.1
    
    for ( q in 1:Q )
      if ( sum(abs(H[PW==q,])) == 0 )
        f[q] = 0
    g[apply(H!=0,2,sum)==0] = 0

    
  # ----------------------------------------------------------------------------------------  
  # Update g
  # Gradients
    gr_H = ((XXUrs-XA+Xtheta)[,Cmap[C2]]*U[,C2[,2]]) %*% Cmat
    gr_g = apply(H*gr_H*sqf[PW],2,sum)

  # Proximal gradient descent with backtracking
    L = sum((A-rep(theta,each=N)-XUrs)^2)
    repeat
    {
      if ( step_g < 1e-7 )
        break
      
      if ( niter <= inititer )
        ng = pmax(g-step_g*gr_g-step_g*lam,0)
      else
        ng = pmax(g-step_g*gr_g-step_g*lam,0)

      G_g = g - ng

      ngH = H*rep(ng,each=P)
      nU = ngH*sqf[PW]
      
      ngHrs = ngH[,C[,1]]*ngH[,C[,2]]
      nUrs = ngHrs*f[PW]
      nXUrs = X%*%nUrs
      nL = sum((A-rep(theta,each=N)-nXUrs)^2)
      
      term2 = sum(gr_g*G_g)
      term3 = sum(G_g*G_g)
      if ( nL < L - term2 + term3/step_g/2 )
      {
        g = ng
        gH = ngH
        U = nU
        gHrs = ngHrs
        Urs = nUrs
        XUrs = nXUrs
        XXUrs = tX%*%XUrs
        break
      }
      step_g = step_g/2
    }
    step_g = step_g*1.1
    
    H[,g==0] = 0
    
  # -------------------------------------------------------------------------------------------  
  # Update f
  # Gradients
    gr_f = rep(0,Q)
    tmp = (XXUrs - XA + Xtheta)*gHrs
    for ( q in 1:Q )
      gr_f[q] = sum(tmp[(PWidx[q]+1):PWidx[q+1],])

  # Proximal gradient descent with backtracking
    L = sum((A-rep(theta,each=N)-XUrs)^2)
    repeat
    {
      if ( step_f < 1e-7 )
        break
      
      if ( niter <= inititer )
        nf = pmax(f-step_f*gr_f-step_f*lam_f/100000,0)
      else
        nf = pmax(f-step_f*gr_f-step_f*lam_f*min(1,(niter-inititer)/gradual),0)

      G_f = f - nf

      sqnf = sqrt(nf)
      nU = gH*sqnf[PW]
      
      nUrs = gHrs*nf[PW]
      nXUrs = X%*%nUrs
      nL = sum((A-rep(theta,each=N)-nXUrs)^2)
      
      term2 = sum(gr_f*G_f)
      term3 = sum(G_f*G_f)
      if ( nL < L - term2 + term3/step_f/2 )
      {
        f = nf
        sqf = sqnf
        U = nU
        Urs = nUrs
        XUrs = nXUrs
        XXUrs = tX%*%XUrs
        break
      }
      step_f = step_f/2
    }
    step_f = step_f*1.1
    
    H[f[PW]==0,] = 0

  # Garbage collection
    gc()
    
  # Convergence check
    if ( niter > inititer+gradual & max(abs(pU-U)) < eps )
      break
    
    print(niter)
  }                                          #end of repeat() 
  
  list(f=f,g=g,H=H,U=U,alpha=alpha,Eomega=Eomega,niter=niter)
}

