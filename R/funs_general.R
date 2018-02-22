

qgmm.marginal2 <- function (u, theta, res = 1000, spread = 5) {
  d <- dim(u)[2]
  m <- length(c(theta$pie))
  n.samples <- round(res * theta$pie)
  n.samples[n.samples == 0] <- 2
  
  # Create grid of evalutation
  s <- NULL
  for (i in 1:d) {
    for (j in 1:m) {
      m.ij <- theta$mu[[j]][i]
      sd.ij <- theta$sigma[[j]][i]
      s <- c(s, seq(m.ij-spread*sd.ij, m.ij+spread*sd.ij, l = n.samples[j]))
    }
  }
  dim(s) <- c(sum(n.samples), d)
  
  # Evaluate on cdf on the grid
  eval <- array(0,dim(s))
  for (i in 1:d){
    for (j in 1:theta$m){
      eval[,i] <- eval[,i] + theta$pie[j]*pnorm(s[,i],mean=theta$mu[[j]][i],sd=theta$sigma[[j]][i])
    }
  }
  
  # Invert function
  z.out <- NULL
  for (j in 1:d) {
    z.out <- c(z.out, approxfun(eval[, j], s[, j], rule = 2)(u[, j]))
  }
  z.out.is.na <- is.na(z.out)
  if (any(z.out.is.na)) {
    z.out[z.out.is.na & u >= 1] <- Inf
    z.out[z.out.is.na & u <= 0] <- -Inf
  }
  dim(z.out) <- c(nrow(u),d)
  return(z.out)
}


getPseduMix <- function(u,para){
  
  Z <- matrix(0,nr=nrow(u), nc=ncol(u)) 
  theta <- vector("list",2)
  
  for (i in 1:2){
    theta[[i]]$m <- 3
    theta[[i]]$d <- 2
    theta[[i]]$pie <- c((1-para[1]),para[1]*c(1-para[i+1],para[i+1]))
    theta[[i]]$mu[[1]] <- rep(0,2)
    theta[[i]]$mu[[2]] <- rep(para[i+5],2)
    theta[[i]]$mu[[3]] <- rep(para[i+3],2)
    theta[[i]]$sigma[[1]] <- rep(1,2)
    theta[[i]]$sigma[[2]] <- rep(1,2)
    theta[[i]]$sigma[[3]] <- rep(para[i+7],2)
    
    Z[,(i*2-1):(i*2)] <- qgmm.marginal2(u[,(i*2-1):(i*2)], theta = theta[[i]],res = 1000,spread = 5)
  }
  return(Z)
}


par_list2vec <- function(par_l){
  par_v <- c(par_l$pi_g[2],par_l$pi_k[,2],par_l$mu,par_l$mu_k0,par_l$sgm,par_l$rho)
  return(par_v)
}


par_vec2list <- function(par_v){
  par_l <- list()
  par_l$pi_g <- c(1-par_v[1],par_v[1])
  par_l$pi_k <- matrix(c(1-par_v[2],1-par_v[3],par_v[2],par_v[3]),ncol=2)
  par_l$mu <- par_v[4:5]
  par_l$mu_k0 <- par_v[6:7]
  par_l$sgm <- par_v[8:9]
  par_l$rho <- par_v[10:11]
  return(par_l)
}

getNegLoglik <- function(par_v,Z){
  
  par <- par_vec2list(par_v)
  m <- length(par$mu)
  
  temp_lik_g1 <- matrix(0, nc=m, nr=nrow(Z))
  for (mi in 1:m){
    temp_lik_g1[,mi] <- par$pi_k[mi,1]*dmvnorm(Z[,c(mi*2-1,mi*2)],mu = rep(par$mu_k0[mi],2),sigma = diag(2)) + 
      par$pi_k[mi,2]*dmvnorm(Z[,c(mi*2-1,mi*2)],mu = rep(par$mu[mi],2),   sigma = sgmToCovm(rep(par$sgm[mi],2),par$rho[mi]))
  }
  temp_lik_g1_1 <- apply(temp_lik_g1,1,prod)
  temp_lik_g0_1 <- dmvnorm(Z,mu = rep(0,4),sigma = diag(4))
  
  neg_log_lik <- -sum(log(par$pi_g[1]*temp_lik_g0_1+par$pi_g[2]*temp_lik_g1_1))
  return(neg_log_lik)
}


sgmToCovm <- function(sgm, rho){
  return(matrix(c(sgm[1]^2,sgm[1]*sgm[2]*rho,sgm[1]*sgm[2]*rho,sgm[2]^2),nc=2,nr=2))
}


Uhat <- function (x) {  # Ranking function
  if (is.vector(x)) {
    x <- matrix(x, length(x), 1)
  }
  apply(x, 2, rank, ties.method = "max")/(nrow(x) + 1)
}


getidrLab <- function(par_v,Z){
  par <- par_vec2list(par_v)
  m <- length(par$mu)
  
  idrm <- matrix(0,nr=dim(Z)[1],nc=2)
  for (mi in 1:m){
    temp0 <- par$pi_k[mi,1]*dmvnorm(Z[,c(mi*2-1,mi*2)],mu = rep(par$mu_k0[mi],2),sigma = diag(2))
    temp1 <- par$pi_k[mi,2]*dmvnorm(Z[,c(mi*2-1,mi*2)],mu = rep(par$mu[mi],2),sigma = sgmToCovm(rep(par$sgm[mi],2),par$rho[mi]))
    idrm[,mi] <- temp0/(temp0+temp1)
  }
  return(idrm)
}

getidrAll <- function(par_v,Z){
  par <- par_vec2list(par_v)
  m <- length(par$mu)
  
  temp10 <- matrix(0, nc=m, nr=nrow(Z))
  for (mi in 1:m){
    temp10[,mi] <- par$pi_k[mi,1]*dmvnorm(Z[,c(mi*2-1,mi*2)],mu = rep(par$mu_k0[mi],2),sigma = diag(2)) + 
      par$pi_k[mi,2]*dmvnorm(Z[,c(mi*2-1,mi*2)],mu = rep(par$mu[mi],2),   sigma = sgmToCovm(rep(par$sgm[mi],2),par$rho[mi]))
  }
  temp1 <- apply(temp10,1,prod)
  temp0 <- dmvnorm(Z,mu = rep(0,4),sigma = diag(4))
  idro <- temp0/(temp1+temp0)
  return(idro)
}















