## Computes the design matrices for LASSO estimation in C
compute_design_C <- function(multi_hawkes,multi_covariates,gamma,beta,alpha,link,observation_matrix) {
  p <- dim(multi_covariates[[1]]$cov[[1]])[1]
  L <- length(multi_covariates[[1]]$times)
  T <- multi_covariates[[1]]$times[L]
  K <- length(multi_hawkes)

  ## Compute Observation Matrix if not provided
  if(is.null(observation_matrix)) {
    Xtilde <- create_observation_matrix(p)
  } else {
    Xtilde <- observation_matrix
  }

  ## Compute Gamma
  Gamma <- matrix(0,nrow=p,ncol=p)
  for(k in 1:K) {
    Gamma <-Gamma+matrix(.Call("compute_gamma",as.integer(p),multi_hawkes[[k]]$EL,as.double(gamma),as.double(T)),ncol=p,nrow=p)
  }
  decompGamma <- eigen(Gamma,symmetric=TRUE)
  if(sum(decompGamma$values<=0)>0) {
    if(sum(decompGamma$values<0)>0) {
      warning("Gamma has negative eigenvalues. In theory this cannot happen. Either there is a mistake in the program or this is due to numerical inaccuracies. This is taken care of in an ad-hoc fashion.")
    }
    ind <- which(decompGamma$values>0)

    eigen_sqrt      <- rep(0,length(decompGamma$values))
    eigen_sqrt[ind] <- sqrt(decompGamma$values[ind])

    eigen_sqrt_inv <- rep(0,length(decompGamma$values))
    eigen_sqrt_inv[ind] <- 1/sqrt(decompGamma$values[ind])
  } else {
    eigen_sqrt     <-   sqrt(decompGamma$values)
    eigen_sqrt_inv <- 1/sqrt(decompGamma$values)
  }
  Gamma_sqrt_inv <- decompGamma$vectors%*%diag(eigen_sqrt_inv)%*%t(decompGamma$vectors)

  ## Compute A
  A <- matrix(0,ncol=p,nrow=p)
  for(k in 1:K) {
    A <- A+matrix(.Call("compute_A",as.integer(p),multi_hawkes[[k]]$EL,as.double(gamma),as.double(T)),ncol=p,nrow=p)
  }

  ## Compute G, V, and v
  G <- matrix(0,ncol=p,nrow=p)
  V <- Matrix::Diagonal(p,x=rep(0,p))
  v <- rep(0,p)

  for(k in 1:K) {
    ## Multiply covariates with beta
    mat <- matrix(.Call("multiply_covariates",multi_covariates[[k]],as.double(beta)),nrow=p)

    ## Apply link function to obtain nu0
    nu0 <- link(mat)

    ## Compute G
    G <- G+matrix(.Call("compute_G",as.integer(p),multi_hawkes[[k]]$EL,as.double(gamma),as.double(T),nu0,multi_covariates[[k]]$times),ncol=p,nrow=p)

    ## Compute Matrix V
    V <- V+Matrix::Diagonal(p,rowSums(t(t(nu0[,-L]^2)*(multi_covariates[[k]]$times[-1]-multi_covariates[[k]]$times[-L]))))

    ## Compute Vector v
    v <- v+.Call("compute_vector_v",multi_hawkes[[k]]$EL,nu0,multi_covariates[[k]]$times)
  }

  ## Compute square root and its inverse of V
  Vsqrt <- Matrix::Diagonal(p,sqrt(Matrix::diag(V)))
  Vsqrt_inv <- Matrix::Diagonal(p,1/sqrt(Matrix::diag(V)))

  ## Compute design and response for LASSO estimation in C
  design_C <- as.matrix(Xtilde%*%decompGamma$vectors%*%diag(eigen_sqrt)%*%t(decompGamma$vectors))
  response_C <- Xtilde%*%Gamma_sqrt_inv%*%(t(A)-t(alpha*G))

  return(list(design=design_C,response=response_C))
}

## Compute coviarance test on a given data set
nodewise_covariance_test <- function(multi_hawkes,multi_covariates,beta,gamma,alpha,estimate_variance=FALSE,link,observation_matrix) {
  ## Get Data
  p <- dim(multi_covariates[[1]]$cov[[1]])[1]

  ## Compute design (and save in variable for simplicity)
  design <- compute_design_C(multi_hawkes,multi_covariates,gamma,beta,alpha,link,observation_matrix)
  X <- design$design
  n <- dim(X)[1]

  ## Compute covariance test for each vertex
  test <- rep(NA,p)
  for(i in 1:p) {
    ## Write response in variable for simplicity
    Y <- design$response[,i]

    ## Compute relevant quantities for test statistics
    U <- t(X)%*%Y
    R <- t(X)%*%X
    m <- max(which(abs(U)>=max(abs(U))))

    ## Compute omega1
    omega1 <- 2*abs(U[m])

    ## Compute omega2
    term1 <- -4*R[m,-m]*(U[m]*U[-m]-R[m,-m]/R[m,m]*U[m]^2)/(R[m,m]*abs(U[m]))
    omega2 <- max((term1-sqrt(term1^2-16*(R[m,-m]^2/R[m,m]^2-1)*(U[-m]-R[m,-m]/R[m,m]*U[m])^2))/(2*R[m,-m]^2/R[m,m]^2-1))

    ## Estimate variance
    if(estimate_variance) {
      out <- glmnet::cv.glmnet(x=X,y=Y)
      estimates <- coef(out)
      residuals <- Y-(estimates[1]+X%*%matrix(estimates[-1],ncol=1))
      sigma <- sd(residuals)
    } else {
      sigma <- sd(Y)
    }

    ## Compute test statistics
    test[i] <- n*omega1*(omega1-omega2)/(4*sigma^2*R[m,m])
  }

  return(test)
}





#### The functions below are typically not directly called by the user.
LASSO_single_line <- function(Y,i,p,T,M_C,omega,m,C.ind.pen) {
  sdY <- sd(Y[,i])*sqrt((m-1)/m)
  sdX <- apply(M_C,2,sd)*sqrt((m-1)/m)
  q <- length(sdX)

  if(sdY==0) {
    ## Y[,i] is identical to zero, in this case the zero vector provides a
    ## perfect solution to the LASSO problem, however, glmnet requires sdY>0
    ## to work properly.
    out <- rep(0,p)

  } else {
    ## Perform LASSO estimation
    K <- q/sum(C.ind.pen/sdX)
    pen.weights <- K*C.ind.pen/sdX

    LASSO <- glmnet::glmnet(t(t(M_C)/sdX),Y[,i]/sdY,intercept=FALSE,standardize=FALSE,lower.limits=rep(0,p),penalty.factor=pen.weights)
    out_raw <- coef(LASSO,s=T*omega[i]/(m*sdY*K),exact=TRUE,x=t(t(M_C)/sdX),y=Y[,i]/sdY,lower.limits=rep(0,p),intercept=FALSE,standardize=FALSE,penalty.factor=pen.weights)[-1]
    out <- sdY*out_raw/sdX
  }

  return(out)
}


compute_lest_squares_theta <- function(par,covariates,C,alpha,hawkes,link) {
  p <- dim(C)[1]
  q <- length(par)-1
  beta <- par[1:q]
  gamma <- par[q+1]
  L <- length(covariates$times)
  T <- covariates$times[L]

  ## Multiply covariates with beta
  mat <- matrix(.Call("multiply_covariates",covariates,as.double(beta)),nrow=p)

  ## Apply link function to obtain nu0
  nu0 <- link(mat)

  ## Compute Matrix V
  V <- Matrix::Diagonal(p,rowSums(t(t(nu0[,-L]^2)*(covariates$times[-1]-covariates$times[-L]))))

  ## Compute Vector v
  v <- .Call("compute_vector_v",hawkes$EL,nu0,covariates$times)

  ## Compute Gamma
  Gamma <- matrix(.Call("compute_gamma",as.integer(p),hawkes$EL,as.double(gamma),as.double(T)),ncol=p,nrow=p)

  ## Compute G
  G <- matrix(.Call("compute_G",as.integer(p),hawkes$EL,as.double(gamma),as.double(T),nu0,covariates$times),ncol=p,nrow=p)

  ## Compute A
  A <- matrix(.Call("compute_A",as.integer(p),hawkes$EL,as.double(gamma),as.double(T)),ncol=p,nrow=p)

  ## Compute Least squares criterion
  LS <- as.numeric(matrix(alpha,nrow=1)%*%V%*%matrix(alpha,ncol=1)+sum(diag(C%*%Gamma%*%t(C)))+2*sum(alpha*diag(C%*%t(G)))-2*sum(alpha*v)-2*sum(diag(C%*%t(A))))

  return(LS)
}

compute_individual_lest_squares_theta <- function(par,covariates,C,alpha,hawkes,link) {
  p <- dim(C)[1]
  q <- length(par)-1
  beta <- par[1:q]
  gamma <- par[q+1]
  L <- length(covariates$times)
  T <- covariates$times[L]

  ## Multiply covariates with beta
  mat <- matrix(.Call("multiply_covariates",covariates,as.double(beta)),nrow=p)

  ## Apply link function to obtain nu0
  nu0 <- link(mat)

  ## Compute Matrix V
  V <- Matrix::Diagonal(p,rowSums(t(t(nu0[,-L]^2)*(covariates$times[-1]-covariates$times[-L]))))

  ## Compute Vector v
  v <- .Call("compute_vector_v",hawkes$EL,nu0,covariates$times)

  ## Compute Gamma
  Gamma <- matrix(.Call("compute_gamma",as.integer(p),hawkes$EL,as.double(gamma),as.double(T)),ncol=p,nrow=p)

  ## Compute G
  G <- matrix(.Call("compute_G",as.integer(p),hawkes$EL,as.double(gamma),as.double(T),nu0,covariates$times),ncol=p,nrow=p)

  ## Compute A
  A <- matrix(.Call("compute_A",as.integer(p),hawkes$EL,as.double(gamma),as.double(T)),ncol=p,nrow=p)

  ## Compute Least squares criterion
  LS <- alpha^2*Matrix::diag(V)+rowSums((C%*%Gamma)*C)+2*alpha*diag(C%*%t(G))-2*alpha*v-2*diag(C%*%t(A))

  return(LS)
}

## This function requires all options from estimate_hawkes other than fit_theta, beta_init, gamma_init.
estimate_hawkes_theta_container <- function(theta,covariates,hawkes,omega,omega_alpha,C.ind.pen,print.level,max_iteration,tol,alpha_init,link,observation_matrix,cluster) {
  ## Read information from data
  p <- dim(C)[1]
  q <- length(theta)-1
  T <- covariates$times[length(covariates$times)]

  ## Compute optimal C and alpha
  opt_theta <- estimate_hawkes(fit_theta=FALSE,beta_init=theta[1:q],gamma_init=theta[q+1],covariates=covariates,hawkes=hawkes,omega=omega,omega_alpha=omega_alpha,C.ind.pen=C.ind.pen,print.level=print.level,max_iteration=max_iteration,tol=tol,alpha_init=alpha_init,link=link,observation_matrix=observation_matrix,cluster=cluster)

  ## Compute objective
  obj <- compute_lest_squares_theta(par=theta,covariates=covariates,C=opt_theta$C,alpha=opt_theta$alpha,hawkes=hawkes,link=link)/T+2*sum(omega*opt_theta$C)+2*omega_alpha*sum(opt_theta$alpha)

  return(obj)
}

estimate_theta_multi_hawkes <- function(theta,multi_covariates,multi_hawkes,omega,omega_alpha,C.ind.pen=NULL,print.level=0,max_iteration=100,tol=0.00001,alpha_init=NULL,link=exp,observation_matrix=NULL,cluster=NULL,return_objective=FALSE) {
  p <- dim(multi_covariates[[1]]$cov[[1]])[1]
  q <- dim(multi_covariates[[1]]$cov[[1]])[2]
  L <- length(multi_covariates[[1]]$times)
  T <- multi_covariates[[1]]$times[L]
  K <- length(multi_hawkes)

  beta <- theta[1:q]
  gamma <- theta[q+1]

  ## Sanity cheks
  if(length(multi_covariates)!=K) {
    stop("'multi_hawkes' and 'multi_covariates' must be of the same length")
  }

  ## Set individual penalties for C estimation to 1 if not provided
  if(is.null(C.ind.pen)) {
    C.ind.pen <- rep(1,p)
  }

  ## Compute Observation Matrix if not provided
  if(is.null(observation_matrix)) {
    Xtilde <- create_observation_matrix(p)
  } else {
    Xtilde <- observation_matrix
  }
  m <- dim(Xtilde)[1]

  ## Set initial values
  if(is.null(alpha_init)) {
    alpha <- rep(1,p)
  } else {
    alpha <- alpha_init
  }
  C <- matrix(NA,ncol=p,nrow=p)

  #### Set dopar if in parallel mode
  if(!is.null(cluster)) {
    `%dopar%` <- foreach::`%dopar%`
  }

  #### Compute design matrices for Lasso estimation
  ## Compute Gamma
  Gamma <- matrix(0,nrow=p,ncol=p)
  for(k in 1:K) {
    Gamma <-Gamma+matrix(.Call("compute_gamma",as.integer(p),multi_hawkes[[k]]$EL,as.double(gamma),as.double(T)),ncol=p,nrow=p)
  }
  decompGamma <- eigen(Gamma,symmetric=TRUE)
  if(sum(decompGamma$values<=0)>0) {
    if(sum(decompGamma$values<0)>0) {
      warning("Gamma has negative eigenvalues. In theory this cannot happen. Either there is a mistake in the program or this is due to numerical inaccuracies. This is taken care of in an ad-hoc fashion.")
    }
    ind <- which(decompGamma$values>0)

    eigen_sqrt      <- rep(0,length(decompGamma$values))
    eigen_sqrt[ind] <- sqrt(decompGamma$values[ind])

    eigen_sqrt_inv <- rep(0,length(decompGamma$values))
    eigen_sqrt_inv[ind] <- 1/sqrt(decompGamma$values[ind])
  } else {
    eigen_sqrt     <-   sqrt(decompGamma$values)
    eigen_sqrt_inv <- 1/sqrt(decompGamma$values)
  }
  Gamma_sqrt_inv <- decompGamma$vectors%*%diag(eigen_sqrt_inv)%*%t(decompGamma$vectors)

  ## Compute design for LASSO estimation in C
  M_C <- as.matrix(Xtilde%*%decompGamma$vectors%*%diag(eigen_sqrt)%*%t(decompGamma$vectors))

  ## Compute A
  A <- matrix(0,ncol=p,nrow=p)
  for(k in 1:K) {
    A <- A+matrix(.Call("compute_A",as.integer(p),multi_hawkes[[k]]$EL,as.double(gamma),as.double(T)),ncol=p,nrow=p)
  }

  ## Compute G, V, and v
  G <- matrix(0,ncol=p,nrow=p)
  V <- Matrix::Diagonal(p,x=rep(0,p))
  v <- rep(0,p)

  for(k in 1:K) {
    ## Multiply covariates with beta
    mat <- matrix(.Call("multiply_covariates",multi_covariates[[k]],as.double(beta)),nrow=p)

    ## Apply link function to obtain nu0
    nu0 <- link(mat)

    ## Compute G
    G <- G+matrix(.Call("compute_G",as.integer(p),multi_hawkes[[k]]$EL,as.double(gamma),as.double(T),nu0,multi_covariates[[k]]$times),ncol=p,nrow=p)

    ## Compute Matrix V
    V <- V+Matrix::Diagonal(p,rowSums(t(t(nu0[,-L]^2)*(multi_covariates[[k]]$times[-1]-multi_covariates[[k]]$times[-L]))))

    ## Compute Vector v
    v <- v+.Call("compute_vector_v",multi_hawkes[[k]]$EL,nu0,multi_covariates[[k]]$times)
  }

  ## Compute square root and its inverse of V
  Vsqrt <- Matrix::Diagonal(p,sqrt(Matrix::diag(V)))
  Vsqrt_inv <- Matrix::Diagonal(p,1/sqrt(Matrix::diag(V)))

  ## Compute Design for LASSO estimation in alpha
  M_alpha <- as.matrix(Xtilde%*%sqrt(V))

  #### Perform Iterative Estimation
  C_old <- C
  alpha_old <- alpha
  par_change <- 0
  iteration <- 1

  TERMINATION_FLAG <- 0
  while(TERMINATION_FLAG==0) {
    #### Estimate C
    ## Compute response for LASSO estimation
    Y <- Xtilde%*%Gamma_sqrt_inv%*%(t(A)-t(alpha*G))

    ## Perform LASSO estimation for each vertex
    if(is.null(cluster)) {
      ## No parallel computation
      for(i in 1:p) {
        C[i,] <- LASSO_single_line(Y,i,p,T,M_C,omega,m,C.ind.pen)
      }
    } else {
      ## Do parallel computations in the provided cluster
      par_out <- foreach::foreach(i=1:p,.combine=rbind,.packages=c('glmnet'),.inorder=FALSE) %dopar% {
        c(i,LASSO_single_line(Y,i,p,T,M_C,omega,m,C.ind.pen))
      }
      ## Bring output in correct order
      C <- par_out[order(par_out[,1]),-1]
    }

    #### Estimate alpha
    ## Compute response for LASSO estimation
    Y <- as.numeric(Xtilde%*%Vsqrt_inv%*%matrix(v-diag(C%*%t(G))))

    ## Perform LASSO estimation for alpha
    sdY <- sd(Y)*sqrt((p-1)/p)
    LASSO <- glmnet::glmnet(M_alpha/sdY,Y/sdY,intercept=FALSE,standardize=FALSE,lower.limits=rep(0,p))
    alpha <- coef(LASSO,s=omega_alpha*p*T/(m*sdY^2),exact=TRUE,x=M_alpha/sdY,y=Y/sdY,lower.limits=rep(0,p),intercept=FALSE,standardize=FALSE)[-1]


    #### Compute progress in alpha and C
    if(iteration==1) {
      C_change <- 2*tol
    } else {
      C_change <- max(abs(C-C_old))
    }
    alpha_change <- max(abs(alpha-alpha_old))
    alphaC_change <- max(c(C_change,alpha_change))

    #### Compute overall progress
    max_change <- max(c(C_change,alpha_change))

    #### Print Status
    if(print.level>0 & iteration>1) {
      cat(sprintf("Finish Iteration %d/%d, c_change=%f, alpha_change=%f Compare to: %f\n",iteration,max_iteration,C_change,alpha_change,tol))
    }
    if(print.level>1) {
      cat("Estimate for C:\n")
      print(C)
      cat("Estimate for alpha:\n")
      print(alpha)
    }

    #### Compute if Termination criterion met
    if(max_change<tol) {
      TERMINATION_FLAG <- 1
    } else if(iteration>=max_iteration) {
      TERMINATION_FLAG <- 1
    }

    C_old <- C
    alpha_old <- alpha
    iteration <- iteration+1
  }

  if(return_objective) {
    return(as.numeric((t(alpha)%*%V%*%alpha+sum(diag(C%*%Gamma%*%t(C)))+2*sum(alpha*diag(C%*%t(G)))-2*sum(alpha*v)-2*sum(diag(C%*%t(A))))/(K*T))+2*sum(omega*C)+2*omega_alpha*sum(alpha))
  } else {
    return(list(C=C,alpha=alpha,beta=beta,gamma=gamma))
  }
}
