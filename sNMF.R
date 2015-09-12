require(matrixStats)

logsumexp <- function(log_p) {
  vec <- log_p - min(log_p)
  return(vec - logSumExp(vec))
}

update.g <- function(g, X, p, c_0, D, Ct, rho) {
  V.seq <- 1:nrow(X)
  K.seq <- 1:dim(p)[1]
  xp <- apply(expand.grid(v=V.seq, k=K.seq), 1, 
              function(itr)sum(X[itr[1], ] * p[itr[2], itr[1], ]))
  xp <- matrix(xp, nrow(X), dim(p)[1])
  return( (1-rho) * g + rho * (c_0/nrow(X) + D/Ct * xp) )
}

update.h <- function(h, a, b, c_0, D, Ct, rho) {
  return( (1-rho) * h + rho * (c_0 + D/Ct * apply(a/b, 1, sum)) )
}

update.a <- function(X, p, a_0) {
  K.seq <- 1:dim(p)[1]
  D.seq <- 1:ncol(X)

  xp <- matrix(
    apply(expand.grid(k=K.seq, d=D.seq), 1, 
          function(itr)sum(X[, itr[2]] * p[itr[1], , itr[2]])), 
    dim(p)[1], ncol(X))
  return(a_0 + xp)
}

update.b <- function(g, h, b_0) {
  return(b_0 + apply(g/h, 2, sum))
}

calc.p <- function(params) {
  E_ln_beta <- digamma(params$g) - log(params$h)
  E_ln_theta <- digamma(params$a) - log(params$b)
  
  V.seq <- 1:nrow(E_ln_beta)
  D.seq <- 1:ncol(E_ln_theta)
  
  p <- apply(expand.grid(v=V.seq, d=D.seq), 1, 
             function(itr)exp(logsumexp(E_ln_beta[itr[1], ] + E_ln_theta[ ,itr[2]])))
  p <- array(p, dim=c(ncol(E_ln_beta),nrow(E_ln_beta),ncol(E_ln_theta)))
  
  return(p)
}

snmf <- function(X, K, a_0=NULL, b_0=NULL, c_0=NULL, 
                 len.Ct=10, t_0=1, kappa=0.5, rho=1, 
                 n.itr=1e2, conv.l=1e2, conv.g=1e2) {
  
  # init variational distribution
  g <- matrix(runif(nrow(X)*K), nrow(X), K)
  h <- runif(K)
  a <- matrix(runif(K*ncol(X)), K, ncol(X))
  b <- runif(K)
  
  # set parameters
  D.seq <- 1:ncol(X)
  params <- list(g=g, h=h, a=a, b=b)
  if (is.null(a_0)) a_0 <- 1.0 / K
  if (is.null(b_0)) b_0 <- 1.0 / K
  if (is.null(c_0)) c_0 <- 0.05 * nrow(X)
  hp <- list(a_0=a_0, b_0=b_0, c_0=c_0)
  
  lb <- vector()  
  for (itr in 1:n.itr) {
    
    #selecting subset of the data
    Ct <- sample(D.seq, len.Ct, replace=FALSE)
    X.sub <- X[, Ct]
    local.params <- list(a=params$a[, Ct], b=params$b)
    
    for (i in 1:1e3) {
      # compute p
      p <- calc.p(list(g=params$g, h=params$h, a=local.params$a, b=local.params$b))
      
      # update local variables
      local.params$a <- update.a(X.sub, p, hp$a_0)
      local.params$b <- update.b(params$g, params$h, hp$b_0)
      
      if(i > 1 && sum(abs(p-pre.p)) < conv.l) break
      pre.p <- p
    }
    
    # update global variables
    params$g <- update.g(params$g, X.sub, p, hp$c_0, ncol(X), len.Ct, rho)
    params$h <- update.h(params$h, local.params$a, local.params$b, hp$c_0, ncol(X), len.Ct, rho)
    params$a[, Ct] <- local.params$a
    params$b <- local.params$b
    
    # update rho
    rho <- (t_0 + itr) ** (-kappa)
    
    lb <- c(lb, sum(abs(X - (params$g/params$h) %*% (params$a/params$b))))
    print(paste(itr, "th iteration : ", lb[itr], sep=""))
    if ( (itr>1) && (abs(lb[itr]-lb[itr-1])<conv.g) ) break
    
  }
  
  plot(lb, type='l')
  
  return(list(params=params, lb=lb))
}

