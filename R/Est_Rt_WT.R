#' Estimate Time-varying Reproduction Numbers with Wallinga-Teunis Method
#'
#'
#'
est.Rt <- function(t = NULL, I_local = NULL, I_import=NULL, 
                   pdf_gi = NULL, h = 5, nsim=100){
  
  mov.sum <- function(x){
    as.vector(stats::filter(x = x, filter = rep(1, h), sides = 2))
  }
  
  est.Rt.sub <- function(t = NULL, I_vect = NULL, Ilc_vect = NULL,
                         pdf_gi = NULL){
    dT <- outer(t, t, FUN = "-")
    dT_lik <- apply(dT, 2, pdf_gi)
    Z_row <- as.vector(dT_lik %*% matrix(I_vect, nrow = Nt, ncol = 1))
    
    Pij <- sapply(1:NCOL(dT_lik), function(j){
      liks_j <- dT_lik[,j]
      (liks_j*Ilc_vect)/Z_row
    })
    
    est.msum <- mov.sum(colSums(Pij, na.rm = T) * I_vect)
    
    I_vect.msum <- mov.sum(I_vect)
    
    est <- est.msum/I_vect.msum
    est[is.na(est) |is.nan(est) | !is.finite(est)] <- 0
    return(list(Pij = Pij, Rt = est))
  }
  
  t <- as.vector(t)
  Nt <- length(t)
  if(is.null(I_import)) I_import <- rep(0, Nt)
  I <- I_local + I_import
  
  # est Rt
  res <- est.Rt.sub(t = t,  I_vect = I, Ilc_vect = I_local, pdf_gi = pdf_gi)
  Rt <- res$Rt
  Pij <- res$Pij


  if(T){
    
    draw.transtree.Rt <- function(P_mat, I_vect, Ilc_vect){
      
      P_mat[is.nan(P_mat)] <- 0
      N0 <- NCOL(P_mat)
      rowsum_P <- rowSums(P_mat, na.rm = T)
      
      sam.tree<- t(sapply(1:NROW(P_mat), function(i){
        if(Ilc_vect[i]==0 | rowsum_P[i]==0) return(rep(0, N0))
        rmultinom(1, Ilc_vect[i], (P_mat[i,]*I_vect))[,1]
      }))
      
      sim.infect.msum <- mov.sum(colSums(sam.tree))
      I_vect.msum <- mov.sum(I_vect)
      
      tmp <- sim.infect.msum/I_vect.msum
      tmp[is.na(tmp)|is.nan(tmp)|!is.finite(tmp)] <- 0
      return(tmp)
    }
    
    sams.Rt <- t(sapply(1:nsim, FUN = function(k){
      draw.transtree.Rt(P_mat = Pij, I_vect = I, Ilc_vect = I_local)
    }))
    Rt_l <- apply(sams.Rt,2,quantile,probs=0.025)
    Rt_u <- apply(sams.Rt,2,quantile,probs=0.975)
    
  }
  
  out <- list(t = t, Rt = Rt, Rt_l = Rt_l, Rt_u = Rt_u)
  
  G <- ggplot(as.data.frame(out))+
    geom_ribbon(aes(x=t, ymin=Rt_l, ymax=Rt_u), alpha=.4)+
    geom_path(aes(x=t, y=Rt))+
    coord_cartesian()
  plot(G)
  
  return(out)
}
