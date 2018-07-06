#' idr model for hierarchical structured data
#'
#' The function measures reproducibility of hierarchical structured data (e.g ChIP-seq data from different labs with multiple replicates).
#'
#' @param expr1 A n*p matrix or data frame of normalized gene expression data. The rows correspond samples (n), the columns correspond genes (p). 
#' 
#' @param expr2 The second gene expression data, should be in the same format, size as expr1.
#' 
#' @param s The tuning parameter for matrixes differences estimation, leave it as NULL to automatically select. 
#' 
#' @param s.seq The candidates for s selection.
#' 
#' @param alpha Prespecified level of false discovery rate. A relatively loose criterion is suggested for determines screening matrix.
#' 
#' @param verbose Set verbose to TURE to show details of s selection.
#' 
#' @details Please refer \bold{``A statistical framework for measuring replicability of heterogeneous omics data from multiple labs'' by Ranalli, Lyu, and Li
}
#' @export


hidr <- function(par0_list,X,
                          nlm_Bound=list(low = c(0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.5, 0.05, 0.05, -0.99, -0.99),
                                         up = c(0.95, 0.95, 0.95,  20,  20,  20,  20,   10,   10,  0.99,  0.99)),
                          nlm_control=list(eval.max=200,iter.max=100,trace=0,rel.tol=0.01),
                          out_control=list(verbose = list(basic=F,par=F),
                                           iterMax = 100,
                                           eps_loglik = 1,
                                           eps_parVec = rep(0.01,length(par0_vec)))){
  u <- Uhat(X)
  par0_vec <- par_list2vec(par0_list)
  epar  <- par0_vec
  epar0 <- rep(0,length(epar))
  
  flag_exit <- F
  itern <- 0 
  neg_loglik_trace <- rep(0,out_control$iterMax) 
  neg_loglik_trace_nlm <- rep(0,out_control$iterMax) 
  
  while(1){
    
    ## updaing psedu data
    
    itern <- itern + 1
    Z <- getPseduMix(u,epar)
    neg_loglik_trace[itern] <- getNegLoglik(epar,Z)
    
    ## verbose displaying
    
    if (out_control$verbose$basic) print(paste(itern,neg_loglik_trace[itern]))
    if (out_control$verbose$par) print(epar)
    
    ## stop conditions
    
    if (itern>=out_control$iterMax) flag_exit <- T
    if (sum(abs(epar-epar0)>out_control$eps_parVec)==0) flag_exit <- T
    if ((itern>1)&&(abs(neg_loglik_trace[itern-1]-neg_loglik_trace[itern])<out_control$eps_loglik)) flag_exit <- T
    if (flag_exit) break
    
    ## updating parameter
    
    output <-  nlminb(epar, getNegLoglik, Z=Z,
                      lower = nlm_Bound$low, upper = nlm_Bound$up,
                      control = nlm_control)
    
    neg_loglik_trace_nlm[itern] <- output$objective 
    epar0 <- epar
    epar <- output$par
  }
  
  idr_lab <- getidrLab(epar,Z)
  idr_all <- getidrAll(epar,Z)
  
  return(list(para=epar, n_inter=itern,  idr_all=idr_all, idr_lab=idr_lab,
              neg_loglik_trace=neg_loglik_trace[which(neg_loglik_trace!=0)],
              neg_loglik_trace_nlm=neg_loglik_trace_nlm[which(neg_loglik_trace!=0)]))
}
