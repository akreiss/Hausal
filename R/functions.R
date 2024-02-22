compute_baseline_intensities <- function(covariates,beta0,alpha,link=exp) {
  p <- dim(covariates$cov[[1]])[1]
  ## Multiply covarites with beta
  mat <- matrix(.Call("multiply_covariates",covariates,as.double(beta0)),nrow=p)

  ## Apply link function
  mat <- link(mat)*alpha

  return(mat)
}

simulate_hawkes <- function(covariates,beta0,gamma,alpha,C,T,link=exp,print.level=1,rand_seed=NA) {
  p <- dim(C)[1]

  ## Compute Baseline
  baseline <- compute_baseline_intensities(covariates,beta0,alpha,link=link)

  #### Step 1: Simulate Baseline Events
  ## Compute Events in a dominating process
  M <- apply(baseline,1,max)
  event_numbers <- rpois(p,T*M)
  proposed_event_times <- runif(sum(event_numbers),min=0,max=T)
  rejection_randomness <- runif(sum(event_numbers),min=0,max=1)

  if(print.level>0) {
    cat("Propose ",sum(event_numbers)," baseline events and perform rejection sampling.\n")
  }

  event_list <- matrix(.Call("baseline_reject",proposed_event_times,as.integer(sum(event_numbers)),as.integer(event_numbers),rejection_randomness,baseline,covariates$times,M),ncol=4)

  #### Step 2: Generate Spin-Off Events
  if(is.na(rand_seed)) {
    rand_seed <- as.integer(Sys.time())
  }
  first_id <- max(event_list[,1])+1
  event_list <- matrix(.Call("generate_spin_off_events",as.double(event_list),as.double(C),rand_seed,as.double(T),as.double(gamma),as.integer(first_id)),ncol=4)

  #### Step 3: Sort according to time
  if(print.level>0) {
    cat("Order Proposed events according to time.\n")
  }
  ord <- order(event_list[,4])
  event_list <- event_list[ord,]

  #### Step 4: Compute intensities
  intint <- .Call("compute_intensity_integrals",event_list,as.double(gamma),covariates$times,baseline,as.double(C))

  ## Prepare Output
  colnames(event_list) <- c("Id","Parent Id","Process","Time")
  out <- list(EL=event_list,intensities=matrix(intint[[1]],nrow=p),integrals=matrix(intint[[2]],nrow=p),baseline=baseline,covariates=covariates)

  return(out)
}

plot_count_intensities <- function(hawkes,T,times) {
  ## Read Information
  p <- dim(hawkes$intensities)[1]
  M <- max(hawkes$intensities)

  ## Plot
  for(i in 1:p) {
    ## Extract events of process i
    ind <- which(hawkes$EL[,3]==i)

    ## Create Plot
    Y <- length(ind)+1
    plot(0,0,type="n",xlim=c(0,T),ylim=c(0,Y),axes=FALSE,main=sprintf("Process %d",i),xlab="Time",ylab="Events")
    axis(1)
    axis(2)

    ## Plot Count process
    if(length(ind)>0) {
      for(k in 1:length(ind)) {
        if(k==1) {
          lines(c(0,hawkes$EL[ind[k],4]),c(0,0))
          if(k==length(ind)) {
            lines(c(hawkes$EL[ind[k],4],T),c(k,k))
          }
        } else if(k==length(ind) & k>1) {
          lines(c(hawkes$EL[ind[k-1],4],hawkes$EL[ind[k],4]),c(k-1,k-1))
          lines(c(hawkes$EL[ind[k],4],T),c(k,k))
        } else {
          lines(c(hawkes$EL[ind[k-1],4],hawkes$EL[ind[k],4]),c(k-1,k-1))
        }

      }
    } else {
      lines(c(0,T),c(0,0))
    }

    ## Plot Intensity
    lines(hawkes$EL[,4],hawkes$intensities[i,]/M*Y,lty=2)
    axis(4,at=0:Y,labels=(0:Y)/Y*M)
  }
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
  V <- Diagonal(p,rowSums(t(t(nu0[,-L]^2)*(covariates$times[-1]-covariates$times[-L]))))

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

estimate_hawkes <- function(covariates,hawkes,omega,omega_alpha,lb,ub,fit_theta=TRUE,print.level=3,max_iteration=100,tol=0.00001,beta_init=NULL,gamma_init=NULL,alpha_init=NULL,link=exp,observation_matrix=NULL) {
  p <- dim(covariates$cov[[1]])[1]
  q <- dim(covariates$cov[[1]])[2]
  L <- length(covariates$times)
  T <- covariates$times[L]

  optimization_args <- list(algorithm="NLOPT_GN_DIRECT_L",xtol_rel=0.0001,print_level=0)

  ## Compute Observation Matrix if not provided
  if(is.null(observation_matrix)) {
    Xtilde <- create_observation_matrix(p)
  } else {
    Xtilde <- observation_matrix
  }
  m <- dim(Xtilde)[1]

  ## Set initial values
  if(is.null(beta_init)) {
    beta <- 0.5*(ub[1:q]+lb[1:q])
  } else {
    beta <- beta_init
  }
  if(is.null(gamma_init)) {
    gamma <- 0.5*(lb[q+1]+ub[q+1])
  } else {
    gamma <- gamma_init
  }
  if(!fit_theta & (is.null(beta_init) | is.null(gamma_init))) {
    stop("If you do not want to fit (beta,gamma), you have to provide initial values for both of them.")
  }
  if(is.null(alpha_init)) {
    alpha <- rep(1,p)
  } else {
    alpha <- alpha_init
  }
  C <- matrix(NA,ncol=p,nrow=p)

  #### Perform Iterative Estimation
  C_old <- C
  alpha_old <- alpha
  par_old <- c(beta,gamma)
  par_change <- 0
  iteration <- 1

  TERMINATION_FLAG <- 0
  BETA_GAMMA_CHANGE <- TRUE
  while(TERMINATION_FLAG==0) {
    #### Compute Auxiliary Matrices and vectors if there was a change in beta and gamma
    if(BETA_GAMMA_CHANGE) {
      ## Compute Gamma
      #      print("Compute Gamma")
      Gamma <- matrix(.Call("compute_gamma",as.integer(p),hawkes$EL,as.double(gamma),as.double(T)),ncol=p,nrow=p)
      decompGamma <- eigen(Gamma,symmetric=TRUE)
      if(sum(decompGamma$values<=0)>0) {
        if(sum(decompGamma$values<0)>0) {
          warning("Gamma has negative eigenvalues. In theory this cannot happen. Either there is a mistake in the program or this is due to numerical inaccuracies. This taken care of in an ad-hoc fashion.")
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
      M_C <- as.matrix(Xtilde%*%diag(eigen_sqrt)%*%t(decompGamma$vectors))

      ## Compute A
      #      print("Compute A")
      A <- matrix(.Call("compute_A",as.integer(p),hawkes$EL,as.double(gamma),as.double(T)),ncol=p,nrow=p)

      ## Multiply covariates with beta
      mat <- matrix(.Call("multiply_covariates",covariates,as.double(beta)),nrow=p)

      ## Apply link function to obtain nu0
      #      print("Compute nu0")
      nu0 <- link(mat)

      ## Compute G
      #      print("Compute G")
      G <- matrix(.Call("compute_G",as.integer(p),hawkes$EL,as.double(gamma),as.double(T),nu0,covariates$times),ncol=p,nrow=p)

      ## Compute Matrix V
      #      print("Compute V")
      V <- Diagonal(p,rowSums(t(t(nu0[,-L]^2)*(covariates$times[-1]-covariates$times[-L]))))
      Vsqrt <- Diagonal(p,sqrt(diag(V)))
      Vsqrt_inv <- Diagonal(p,1/sqrt(diag(V)))

      ## Compute Vector v
      #      print("Compute vector v")
      v <- .Call("compute_vector_v",hawkes$EL,nu0,covariates$times)

      ## Compute Design for LASSO estimation in alpha
      M_alpha <- as.matrix(Xtilde%*%sqrt(V))

      ## Recompute the above only if there was a change in beta and gamma
      BETA_GAMMA_CHANGE <- FALSE
    }


    #### Estimate C
    ## Compute response for LASSO estimation
    Y <- Xtilde%*%t(decompGamma$vectors)%*%Gamma_sqrt_inv%*%(t(A)-t(alpha*G))

    ## Perform LASSO estimation for each vertex
    for(i in 1:p) {
      sdY <- sd(Y[,i])*sqrt((p-1)/p)
      LASSO <- glmnet(M_C/sdY,Y[,i]/sdY,intercept=FALSE,standardize=FALSE,lower.limits=rep(0,p))
      C[i,] <- coef(LASSO,s=omega[i]*p*T/(m*sdY^2))[-1]
    }


    #### Estimate alpha
    ## Compute response for LASSO estimation
    Y <- as.numeric(Xtilde%*%Vsqrt_inv%*%matrix(v-diag(C%*%t(G))))

    ## Perform LASSO estimation for alpha
    sdY <- sd(Y)*sqrt((p-1)/p)
    LASSO <- glmnet(M_alpha/sdY,Y/sdY,intercept=FALSE,standardize=FALSE,lower.limits=rep(0,p))
    alpha <- coef(LASSO,s=omega_alpha*p*T/(m*sdY^2))[-1]


    ##### Estimate theta=(beta,gamma) if asked to do so
    if(fit_theta) {
      if(iteration %% 10==1) {
        out <- nloptr(c(beta,gamma),compute_lest_squares_theta,opts=optimization_args,ub=ub,lb=lb,covariates=covariates,C=C,alpha=alpha,hawkes=hawkes,link=link)
        beta <- out$solution[1:q]
        gamma <- out$solution[q+1]

        BETA_GAMMA_CHANGE <- TRUE

        par_change <- max(abs(c(beta,gamma)-par_old))
        par_old <- c(beta,gamma)
      }
    }

    #### Compute Progress
    if(iteration==1) {
      C_change <- 2*tol
    } else {
      C_change <- max(abs(C-C_old))
    }
    alpha_change <- max(abs(alpha-alpha_old))
    max_change <- max(c(C_change,alpha_change,par_change))

    #### Print Status
    if(print.level>0 & iteration>1) {
      cat(sprintf("Finish Iteration %d/%d, c_change=%f, alpha_change=%f, par_change=%f Compare to: %f\n",iteration,max_iteration,C_change,alpha_change,par_change,tol))
    }
    if(print.level>1) {
      cat("Estimate for C:\n")
      print(C)
      cat("Estimate for alpha:\n")
      print(alpha)
      cat("Estimate for par:\n")
      print(c(beta,gamma))
    }

    #### Compute if Termination criterion met
    if(max_change<tol) {
      TERMINATION_FLAG <- 1
    } else if(iteration>max_iteration) {
      TERMINATION_FLAG <- 1
    }

    C_old <- C
    alpha_old <- alpha
    iteration <- iteration+1
  }

  return(list(C=C,alpha=alpha,beta=beta,gamma=gamma))
}

debias_Hawkes <- function(covariates,hawkes,est_hawkes,link=exp,observation_matrix=NULL) {
  p <- length(est_hawkes$alpha)
  q <- length(est_hawkes$beta)
  L <- length(covariates$times)
  T <- covariates$times[L]

  ## Check if Design Matrix is provided
  if(is.null(observation_matrix)) {
    tildeX <- create_observation_matrix(1+q+p+p^2-1)
  } else {
    tildeX <- observation_matrix
  }

  #### Compute Auxiliary Information
  ## Multiply covariates with beta
  mat <- matrix(.Call("multiply_covariates",covariates,as.double(est_hawkes$beta)),nrow=p)

  ## Apply link function to obtain nu0
  nu0 <- link(mat)

  ## Compute Matrix V
  V_diag <- rowSums(t(t(nu0[,-L]^2)*(covariates$times[-1]-covariates$times[-L])))

  ## Compute V Matrix derivatives
  Vderivs <- .Call("dbeta_V",nu0,covariates)
  Vderivs[[1]] <- matrix(Vderivs[[1]],ncol=q)
  for(i in 1:q) {
    Vderivs[[2]][[i]] <- matrix(Vderivs[[2]][[i]],ncol=q)
  }

  ## Compute Vector v
  v <- .Call("compute_vector_v",hawkes$EL,nu0,covariates$times)

  ## Compute derivatives of v
  vecvderivs <- .Call("dbeta_vector_v",hawkes$EL,nu0,covariates)
  vecvderivs[[1]] <- matrix(vecvderivs[[1]],ncol=q)
  for(i in 1:q) {
    vecvderivs[[2]][[i]] <- matrix(vecvderivs[[2]][[i]],ncol=q)
  }

  ## Compute G
  G <- matrix(.Call("compute_G",as.integer(p),hawkes$EL,as.double(est_hawkes$gamma),as.double(T),nu0,covariates$times),ncol=p,nrow=p)

  ## Compute derivatives of G
  Gderivs <- .Call("compute_deriv_G",as.integer(p),hawkes$EL,as.double(est_hawkes$gamma),nu0,covariates)
  Gderivs[[1]] <- matrix(Gderivs[[1]],ncol=p)
  Gderivs[[2]] <- matrix(Gderivs[[2]],ncol=p)
  for(i in 1:q) {
    Gderivs[[3]][[i]] <- matrix(Gderivs[[3]][[i]],ncol=p)
    Gderivs[[4]][[i]] <- matrix(Gderivs[[4]][[i]],ncol=p)
    for(j in 1:q) {
      Gderivs[[5]][[i]][[j]] <- matrix(Gderivs[[5]][[i]][[j]],ncol=p)
    }
  }


  ## Compute Gamma
  Gamma <- matrix(.Call("compute_gamma",as.integer(p),hawkes$EL,as.double(est_hawkes$gamma),as.double(T)),ncol=p,nrow=p)

  ## Compute derivatives of Gamma
  Gammaderivs <- .Call("compute_gamma_deriv",as.integer(p),hawkes$EL,as.double(est_hawkes$gamma),as.double(T))
  Gammaderivs[[1]] <- matrix(Gammaderivs[[1]],ncol=p)
  Gammaderivs[[2]] <- matrix(Gammaderivs[[2]],ncol=p)

  ## Compute A
  A <- matrix(.Call("compute_A",as.integer(p),hawkes$EL,as.double(est_hawkes$gamma),as.double(T)),ncol=p,nrow=p)

  ## Compute derivative of A
  Aderivs <- .Call("compute_A_deriv",as.integer(p),hawkes$EL,as.double(est_hawkes$gamma),as.double(T))
  Aderivs[[1]] <- matrix(Aderivs[[1]],ncol=p);
  Aderivs[[2]] <- matrix(Aderivs[[2]],ncol=p);


  #### Compute Gradient
  grad <- .Call("compute_derivatives",as.double(V_diag),as.double(est_hawkes$alpha),as.double(v),as.double(est_hawkes$C),as.double(G),as.double(Gamma),as.double(A),as.integer(q))/(p*T)

  #### Compute Sigma
  Sigma <- matrix(NA,nrow=q+1+p+p^2,ncol=q+1+p+p^2)

  for(k in 1:q) {
    for(l in k:q) {
      ## Second derivatives with respect to beta
      Sigma[k,l] <- (sum(est_hawkes$alpha^2*Vderivs[[2]][[k]][,l])+2*sum(alpha*diag(est_hawkes$C%*%t(Gderivs[[5]][[k]][[l]])))-2*sum(alpha*vecvderivs[[2]][[k]][,l]))/(p*T)
      Sigma[l,k] <- Sigma[k,l]
    }
    ## Derivative with respect to gamma and beta
    Sigma[q+1,k] <- 2*sum(alpha*diag(C%*%t(Gderivs[[4]][[k]])))/(p*T)
    Sigma[k,q+1] <- Sigma[q+1,k]

    ## Derivatives with respect to alpha and beta
    Sigma[k,(q+2):(q+1+p)] <- as.numeric(2*Diagonal(p,Vderivs[[1]][,k])%*%est_hawkes$alpha+2*diag(C%*%t(Gderivs[[3]][[k]]))-2*vecvderivs[[1]][,k])/(p*T)
    Sigma[(q+2):(q+1+p),k] <- Sigma[k,(q+2):(q+1+p)]

    ## Derivatives with respect to C and beta
    Sigma[k,(q+1+p+1):(q+1+p+p^2)] <- as.numeric(2*Diagonal(p,est_hawkes$alpha)%*%Gderivs[[3]][[k]])/(p*T)
    Sigma[(q+1+p+1):(q+1+p+p^2),k] <- Sigma[k,(q+1+p+1):(q+1+p+p^2)]
  }

  ## Second Derivative with respect to gamma
  Sigma[q+1,q+1] <- (sum(diag(est_hawkes$C%*%Gammaderivs[[2]]%*%t(est_hawkes$C)))+2*sum(est_hawkes$alpha*diag(est_hawkes$C%*%t(Gderivs[[2]])))-2*sum(diag(est_hawkes$C%*%t(Aderivs[[2]]))))/(p*T)

  ## Derivative with respect to gamma and alpha
  Sigma[q+1,(q+2):(q+1+p)] <- 2*diag(est_hawkes$C%*%t(Gderivs[[1]]))/(p*T)
  Sigma[(q+2):(q+1+p),q+1] <- Sigma[q+1,(q+2):(q+1+p)]

  ## Derivative with respect to gamma and C
  Sigma[q+1,(q+1+p+1):(q+1+p+p^2)] <- as.numeric(est_hawkes$C%*%Gammaderivs[[1]]+est_hawkes$C%*%t(Gammaderivs[[1]])+2*Diagonal(p,est_hawkes$alpha)%*%Gderivs[[1]]-2*Aderivs[[1]])/(p*T)
  Sigma[(q+1+p+1):(q+1+p+p^2),q+1] <- Sigma[q+1,(q+1+p+1):(q+1+p+p^2)]

  ## Second derivative with respect to alpha
  Sigma[(q+2):(q+1+p),(q+2):(q+1+p)] <- 2*diag(V_diag)/(p*T)

  ## Compute derivative with respect to alpha and C
  Sigma[(q+2):(q+1+p),(q+1+p+1):(q+1+p+p^2)] <- 2*matrix(.Call("compute_alphaC_deriv",G),nrow=p)/(p*T)
  Sigma[(q+1+p+1):(q+1+p+p^2),(q+2):(q+1+p)] <- t(Sigma[(q+2):(q+1+p),(q+1+p+1):(q+1+p+p^2)])

  ## Compute Second Derivative with respect to C
  Sigma[(q+1+p+1):(q+1+p+p^2),(q+1+p+1):(q+1+p+p^2)] <- matrix(.Call("compute_d2C_deriv",Gamma),nrow=p^2)/(p*T)

  #### Compute Nodewise LASSO for the first columns of Sigma corresponding to beta and gamma
  ## Compute sqrt of Sigma
  Sigma_eigen <- eigen(Sigma*(p*T))
  Lambda <- Diagonal(1+q+p+p^2,Sigma_eigen$values)

  ## Compute Nodewise LASSO
  Theta <- matrix(NA,ncol=1+q+p+p^2,nrow=q+1)
  for(j in 1:(q+1)) {
    X <- t(t(Sigma_eigen$vectors)[,-j])%*%Lambda%*%t(Sigma_eigen$vectors)
    Z <- as.numeric(tildeX%*%X[,j])
    M <- as.matrix(tildeX%*%X[,-j])

    node_wise_lasso <- cv.glmnet(M,Z,penalty.factor=c(rep(0,q),rep(1,p+p^2)),intercept=FALSE,nfolds=5,standardize=FALSE)
    vec <- coef(node_wise_lasso)[-1]

    tau <- as.numeric(matrix(Sigma_eigen$vectors[j,],nrow=1)%*%Lambda%*%(matrix(t(Sigma_eigen$vectors)[,j],ncol=1)-t(Sigma_eigen$vectors)[,-j]%*%vec))/(p*T)

    Theta[j,-j] <- -vec/tau
    Theta[j,j] <- 1/tau
  }

  ## Compute De-Biased Estimator
  theta_debiased <- c(est_hawkes$beta,est_hawkes$gamma)+Theta%*%matrix(grad,ncol=1)

  return(list(grad=grad,Sigma=Sigma,Theta=Theta,beta_debiased=theta_debiased[1:q],gamma_debiased=theta_debiased[q+1]))
}

## Computes the complete estimator with all stages
NetHawkes <- function(covariates,hawkes,omega,omega_alpha,lb,ub,print.level=2,max_iteration=100,tol=0.00001,link=exp,observation_matrix_network=NULL,observation_matrix_debiasing=NULL) {
  ## Perform first stage estimation
  if(print.level>0) {
    cat("Perform the first stage estimation.\n")
  }
  est_first_stage <- estimate_hawkes(covariates=covariates,hawkes=hawkes,omega=omega,omega_alpha=omega_alpha,lb=lb,ub=ub,fit_theta=TRUE,print.level=print.level,max_iteration=max_iteration,tol=tol,beta_init=NULL,gamma_init=NULL,alpha_init=NULL,link=link,observation_matrix=observation_matrix_network)

  ## Debiasing
  if(print.level>0) {
    cat("Debias the first stage estimator.\n")
  }
  debiased_est <- debias_Hawkes(covariates=covariates,hawkes=hawkes,est_hawkes=est_first_stage,link=link,observation_matrix=observation_matrix_debiasing)

  ## Compute Network estimate with debiased estimator
  if(print.level>0) {
    cat("Compute second stage estimator.\n")
  }
  est_second_stage <- estimate_hawkes(covariates=covariates,hawkes=hawkes,omega=omega,omega_alpha=omega_alpha,lb=lb,ub=ub,fit_theta=FALSE,print.level=print.level,max_iteration=max_iteration,tol=tol,beta_init=debiased_est$beta_debiased,gamma_init=debiased_est$gamma_debiased,alpha_init=est_first_stage$alpha,link=link,observation_matrix=observation_matrix_network)

  return(list(first_stage=est_first_stage,second_stage=est_second_stage,debiasing=debiased_est))
}


## Computes the observation matrices Xtilde and Xbar for dimension n
## Input:
## n - The dimension of the required matrix
## Output: List of the following two matrices:
## tildeX    - Matrix \tilde{X}_n of dimension 3n/2 x n or (3n+1)/2 x n.
## tildeXinv - Matrix \overline{X}_n of dimension n x n.
create_observation_matrix <- function(n) {
  if(n %% 2==0) {
    ev <- n/2
    od <- n/2
  } else {
    ev <- (n-1)/2
    od <- (n-1)/2+1
  }
  M <- 2*ev+3*od
  i <- 1:M
  j <- 1:M
  x <- 1:M

  ## Even Diagonal
  ind <- 1:ev
  i[ind] <- 2*(1:ev)
  j[ind] <- 2*(1:ev)
  x[ind] <- 1/sqrt(2)

  ## Odd Diagonal
  ind <- (ev+1):n
  i[ind] <- 2*(1:od)-1
  j[ind] <- 2*(1:od)-1
  x[ind] <- 1/sqrt(6)

  ## Upper Off-diagonal
  ind <- (n+1):(n+od)
  i[ind] <- 2*(1:od)-1
  j[ind] <- 2*(1:od)
  x[ind] <- 1/sqrt(6)

  ## Lower Off-diagonal
  ind <- (n+od+1):(n+od+ev)
  i[ind] <- 2*(1:ev)
  j[ind] <- 2*(1:ev)-1
  x[ind] <- -1/sqrt(2)

  ## Row Corrections
  ind <- (n+od+ev+1):(n+2*od+ev)
  i[ind] <- 2*(1:od)-1
  j[ind] <- n+1:od
  x[ind] <- -2/sqrt(6)

  ## Create Matrix
  M <- sparseMatrix(i=i,j=j,x=x)

  ## Last row requires modification for odd n
  if(n %%2 ==1) {
    M[n,n] <- 1/sqrt(2)
    M[n,n+1] <- 0
    M[n,n+od] <- -1/sqrt(2)
  }

  return(t(M))
}

plot_interactions <- function(estHawkes,vertex.scaling=1,edge.scaling=1,vertex.names=NULL,show.plot=TRUE,...) {
  ## Create Network
  G <- graph_from_adjacency_matrix(edge.scaling*estHawkes$C,mode="directed",weighted="weight")

  ## Change Vertex names if applicable
  if(!is.null(vertex.names)) {
    V(G)$name <- vertex.names
  }

  ## Create Plot if asked
  if(show.plot) {
    ## Generate Node Sizes
    node_sizes <- vertex.scaling*(0.1+sqrt(estHawkes$alpha))

    ## Plot
    plot(G,vertex.size=node_sizes,edge.width=E(G)$weight,...)
  }

  return(G)
}
