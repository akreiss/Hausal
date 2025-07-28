#' Simulate an Event Network from provided covariates
#'
#' Provides a simulated event history from the Hawkes Causal Model. The true
#' model parameters and covariate processes have to be provided. In addition to
#' the simulated events, `simulate_hawkes` also returns the overall individual
#' intensity functions, the baseline intensities and the parts of the intensity
#' due to mutual excitation.
#'
#' The result of the function is random and will therefore be different for
#' different calls unless when the random seed is set before calling the
#' function and setting the `rand_seed` value.
#'
#' `covariates` is a list of two elements `times` and `cov`
#' * `times` is a vector of increasing time points with last element equal to
#'   `T`. `i`-th element of `cov` contains the covariate information valid for
#'   the interval `covariates$times[i]` till `covariates$times[i+1]`.
#' * `cov` is itself a list. Its length equals the length of `times`. Each
#'   element of `cov` is a pxq-matrix, where p is the number of vertices in the
#'   network and q is the number of covariates. Each row of this matrix
#'   corresponds to the covariate vector valid for the respective vertex.
#'
#' Details of the simulation process are provided in our paper.
#'
#' @param covariates Covariate information, its exact format is described in the
#'   Details below.
#' @param beta0 Vector of length q, the true parameter of the baseline
#'   intensity.
#' @param gamma Decay rate (typically positive) in the excitation kernels.
#' @param alpha Non-negative vector of length p (number of vertices), where each
#'   entry gives the individual activity of the corresponding vertex.
#' @param C Weighted adjacency matrix of the network, must be non-negative and
#'   of dimension p x p.
#' @param T Positive number giving the end of observation period.
#' @param link A function that can be applied to vectors returning a
#'   non-negative vector of the same length, it will be applied to linear
#'   transformation of the covariates with `beta0` to obtain the baseline
#'   intensity. The defualt choice is  `exp`.
#' @param print.level A single number, if 0 (the default) no status information
#'   is printed, if positive some updates will be given about the process of the
#'   the simulation.
#' @param rand_seed Integer, if `NA` (the default) the current time stamp is
#'   converted to an integer value. For the simulation it is necessary to
#'   compute random numbers in a C sub-routine, `rand_seed` will be passed to
#'   this routine. See also the Details section.
#'
#' @returns A list with the following elements:
#' * `EL`: A matrix with four columns, each row corresponds to an event. The
#'   third column specfies in which vertex the event happened and the fourth
#'   column at which time. The rows are ordered according to the fourth column.
#'   The first column contains a unique id (starting from 1) which identifies
#'   the event. The second column specifies from which parent event the current
#'   event was generated, events generated through the baseline have a 0 here.
#' * `intensities`: A matrix with p rows and each column corresponds to an event
#'   (in the same order as in `EL`). Each column gives the intensities of the
#'   corresponding processes at the time of the event.
#' * `integrals`: A matrix with p rows and each column corresponds to an event
#'   (in the same order as in `EL`). Each column provides the values of the
#'   integrals over the excitation kernels with respect to the corresponding
#'   counting process (per vertex) and at the time of the corresponding event.
#' * `baseline`: A matrix with p rows and number of columns equal to the length
#'   of `covariates$times`. The i-th column contains the baseline intensities
#'   of the corresponding vertices at time `covariates$times[i]`.
#' * `covariates`: A copy of the input `covariates`.
#'
#' @export
simulate_hawkes <- function(covariates,beta0,gamma,alpha,C,T,link=exp,print.level=0,rand_seed=NA) {
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
  event_list <- event_list[ord,,drop=FALSE]

  #### Step 4: Compute intensities
  intint <- .Call("compute_intensity_integrals",event_list,as.double(gamma),covariates$times,baseline,as.double(C))

  ## Prepare Output
  colnames(event_list) <- c("Id","Parent Id","Process","Time")
  out <- list(EL=event_list,intensities=matrix(intint[[1]],nrow=p),integrals=matrix(intint[[2]],nrow=p),baseline=baseline,covariates=covariates)

  return(out)
}

#' Compute the Baseline Intensity
#'
#' When provided with the covariates, `compute_baseline_intensities` returns the
#' baseline intensities at the times when the covariates change.
#'
#' @inheritParams simulate_hawkes
#' @param covariates Covariate information, its exact format is described in the
#'   Details of [simulate_hawkes()].
#'
#' @returns `compute_baseline_intensities` returns a matrix wit p (number of
#'   vertices) rows and the number of columns equals the length of
#'   `covariates$times`. The i-th column of the output gives the baseline
#'   intensities of all vertices at the time `covariates$times[i]`.
#' @export
compute_baseline_intensities <- function(covariates,beta0,alpha,link=exp) {
  p <- dim(covariates$cov[[1]])[1]
  ## Multiply covarites with beta
  mat <- matrix(.Call("multiply_covariates",covariates,as.double(beta0)),nrow=p)

  ## Apply link function
  mat <- link(mat)*alpha

  return(mat)
}



#' Fits a Network Hawkes Model with Covariates
#'
#' `estimate_hawkes` estimates the parameter of the Hawkes Causal Model as
#' described in the paper. It requires a covariate process, cf.
#' [simulate_hawkes()] for a description, the LASSO parameters and lower and
#' upper bounds on the parameters as inputs. The observed point process is
#' provided as a list in the same format as returned by [simulate_hawkes()]
#'
#' The estimation is an iterative procedure as mentioned in the paper. In each
#' iteration, we first update `C`, then `alpha` and then `beta` and `gamma`.
#' Therefore, initial values are required for all parameters other than `C`, cf.
#' the corresponding parameters of the function. These iterations are repeated
#' until the maximum number of iterations is reached or until the change of the
#' parameters after one iteration lies below the threshold. It is possible to
#' provide values for `beta` and `gamma` and to estimate only `C` and `alpha`.
#'
#' @inheritParams simulate_hawkes
#' @param covariates Covariate information, its exact format is described in the
#'   Details of [simulate_hawkes()].
#' @param hawkes The observed point process. This is a list with at least the
#'   element `EL`. This element needs to have the same matrix form as in the
#'   output of [simulate_hawkes()].
#' @param omega Vector of length p (number of vertices), containing non-negative
#'   tuning parameter for the LASSO penalty for the corresponding row of the
#'   adjacency matrix.
#' @param omega_alpha A non-negative number, containing the tuning parameter for
#'   the LASSO penalty with respect to the individuals intensities a.
#' @param lb,ub Vectors of length q+1, q equals the dimension of covariates. The
#'   first q entries of `lb` and `ub` provide lower and upper bounds on beta,
#'   respectively. The last entry provides a lower (resp. upper) bound on gamma.
#' @param C.ind.pen Individual weights for the LASSO estimation of C, which are
#'   passed to glmnet as penalty.factor. The weights are internally rescaled. To
#'   achieve the criterion function from the paper, set it equal to `NULL` (the
#'   default), then, all vertices receive the same weight.
#' @param fit_theta Logical value, if TRUE (the default) the parameters beta and
#'   gamma are also fitted. If FALSE, beta and gamma are fixed equal to the
#'   provided values in `beta_init`, `gamma_init`.
#' @param print.level Integer which specifies how much information about the
#'   iteration should be printed: 0 (the default) means no information, 1
#'   provides some information and any number larger than 1 results in printing
#'   of all intermediate steps. Note that the nloptr print.level is always equal
#'   to 0.
#' @param max_iteration Maximal number of iterations after which the iteration
#'   is stopped. Default is 100.
#' @param tol When the parameters have changed after one iteration less than the
#'   value provided in `tol`, the iteration is stopped and the current value is
#'   returned as the result. The default is 0.00001.
#' @param beta_init,gamma_init If `fit_theta=FALSE` these parameters provide
#'   starting values for beta and gamma, respectively. If `NULL` (the default),
#'   the middle between `lb` and `ub` is chosen. If `fit_theta=FALSE`, see
#'   `fit_theta`. In this case both `beta_init` and `gamma_init` must be
#'   specified.
#' @param alpha_init Vector of initial values for alpha, unless it is `NULL`
#'   (the default) in which case all alphas are initialized with 1.
#' @param observation_matrix If `NULL` (the default), the observation matrix is
#'   computed within `estimate_hawkes`. For a single call of `estimate_hawkes`
#'   this makes no difference. However, for repeated calls it might be more
#'   efficient to compute the matrix using [create_observation_matrix()] and
#'   pass it as a parameter here.
#' @param cluster If `NULL` (the default) serial computations are executed. If a
#'   cluster as returned by `makeCluster` is provided (after calling
#'   `registerDoParallel(cluster)`), the estimation of C is executed in
#'   parallel. This requires the packages `parallel`, `doParallel`, and
#'   `foreach`.
#'
#' @returns Returns a list with the elements `C`, `alpha`, `beta`, and `gamma`
#'   which contain the estimates for the respective parameters.
#'
#' @export
estimate_hawkes <- function(covariates,hawkes,omega,omega_alpha,lb,ub,C.ind.pen=NULL,fit_theta=TRUE,print.level=0,max_iteration=100,tol=0.00001,beta_init=NULL,gamma_init=NULL,alpha_init=NULL,link=exp,observation_matrix=NULL,cluster=NULL) {
  p <- dim(covariates$cov[[1]])[1]
  q <- dim(covariates$cov[[1]])[2]
  L <- length(covariates$times)
  T <- covariates$times[L]

  optimization_args <- list(algorithm="NLOPT_GN_DIRECT_L",xtol_rel=0.0001,print_level=0)

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

  #### Set dopar if in parallel mode
  if(!is.null(cluster)) {
    `%dopar%` <- foreach::`%dopar%`
  }

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
      Gamma <- matrix(.Call("compute_gamma",as.integer(p),hawkes$EL,as.double(gamma),as.double(T)),ncol=p,nrow=p)
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
      V <- Matrix::Diagonal(p,rowSums(t(t(nu0[,-L]^2)*(covariates$times[-1]-covariates$times[-L]))))
      Vsqrt <- Matrix::Diagonal(p,sqrt(Matrix::diag(V)))
      Vsqrt_inv <- Matrix::Diagonal(p,1/sqrt(Matrix::diag(V)))

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


    ##### Estimate theta=(beta,gamma) if asked to do so
    if(fit_theta) {
      if(alphaC_change<tol | iteration %% 10==1 | iteration==max_iteration) {
        if(print.level>1) {
          cat("Refit theta\n")
        }
        out <- nloptr::nloptr(c(beta,gamma),compute_lest_squares_theta,opts=optimization_args,ub=ub,lb=lb,covariates=covariates,C=C,alpha=alpha,hawkes=hawkes,link=link)
        beta <- out$solution[1:q]
        gamma <- out$solution[q+1]

        BETA_GAMMA_CHANGE <- TRUE

        par_change <- max(abs(c(beta,gamma)-par_old))
        par_old <- c(beta,gamma)
      }
    }

    #### Compute overall progress
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
    } else if(iteration>=max_iteration) {
      TERMINATION_FLAG <- 1
    }

    C_old <- C
    alpha_old <- alpha
    iteration <- iteration+1
  }

  return(list(C=C,alpha=alpha,beta=beta,gamma=gamma))
}

#' De-Bias a Given Estimator
#'
#' `debias_hawkes` takes an estimator, e.g., computed through
#' [estimate_hawkes()] and performs the de-biasing on the estimates for `beta`
#' and `gamma`.
#'
#' De-Biasing is required in order to remove the bias from `beta` and `gamma`
#' which is introduced in [estimate_hawkes()] due to the LASSO penalty.
#'
#' @inheritParams estimate_hawkes
#' @param est_hawkes An estimated Hawkes Causal Model, the estimate has to be
#'   formatted as the output of [estimate_hawkes()].
#' @param debias_thresh Passed to glmnet(), the default value is `1e-14`,
#'   glmnet() suggests `1e-07`. Changing this value can speed up the
#'   computations.
#' @param debias_maxit Passed to glmnet(), the default value is `1e+08`,
#'   glmnet() suggests `1e+05`. Changing this value can speed up the
#'   computations.
#' @param debias_exact If `TRUE`, the default, the node-wise lasso is computed
#'   exactly by coef() from the glmnet package. If `FALSE`, an approximation is
#'   used, which can speed-up computations.
#'
#' @returns `debias_hawkes` returns a list with the following elements:
#'  * `grad`: A vector containing the derivative of the criterion function.
#'  * `Sigma`: A matrix containing the second derivative of the criterion
#'    function.
#'  * `Theta`: Matrix with q+1 rows (q being the dimension of the covariates)
#'    containing the first q+1 rows of the approximation of the inverse of Sigma
#'    via node-wise LASSO.
#'  * `beta_debiased`: De-biased estimate for beta.
#'  * `gamma_debiased`: De-biased estimate for gamma.
#'  * `tuning_parameter`: The value that is used for all sigma_j
#'
#' @export
debias_Hawkes <- function(covariates,hawkes,est_hawkes,link=exp,observation_matrix=NULL,debias_thresh=1e-14,debias_maxit=100000000,debias_exact=TRUE) {
  p <- length(est_hawkes$alpha)
  q <- length(est_hawkes$beta)
  L <- length(covariates$times)
  T <- covariates$times[L]

  ## Check if Design Matrix is provided
  if(is.null(observation_matrix)) {
    tildeX <- create_observation_matrix(1+q+p+p^2)
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
      Sigma[k,l] <- (sum(est_hawkes$alpha^2*Vderivs[[2]][[k]][,l])+2*sum(est_hawkes$alpha*diag(est_hawkes$C%*%t(Gderivs[[5]][[k]][[l]])))-2*sum(est_hawkes$alpha*vecvderivs[[2]][[k]][,l]))/(p*T)
      Sigma[l,k] <- Sigma[k,l]
    }
    ## Derivative with respect to gamma and beta
    Sigma[q+1,k] <- 2*sum(est_hawkes$alpha*diag(est_hawkes$C%*%t(Gderivs[[4]][[k]])))/(p*T)
    Sigma[k,q+1] <- Sigma[q+1,k]

    ## Derivatives with respect to alpha and beta
    Sigma[k,(q+2):(q+1+p)] <- as.numeric(2*Matrix::Diagonal(p,Vderivs[[1]][,k])%*%est_hawkes$alpha+2*diag(est_hawkes$C%*%t(Gderivs[[3]][[k]]))-2*vecvderivs[[1]][,k])/(p*T)
    Sigma[(q+2):(q+1+p),k] <- Sigma[k,(q+2):(q+1+p)]

    ## Derivatives with respect to C and beta
    Sigma[k,(q+1+p+1):(q+1+p+p^2)] <- as.numeric(2*Matrix::Diagonal(p,est_hawkes$alpha)%*%Gderivs[[3]][[k]])/(p*T)
    Sigma[(q+1+p+1):(q+1+p+p^2),k] <- Sigma[k,(q+1+p+1):(q+1+p+p^2)]
  }

  ## Second Derivative with respect to gamma
  Sigma[q+1,q+1] <- (sum(diag(est_hawkes$C%*%Gammaderivs[[2]]%*%t(est_hawkes$C)))+2*sum(est_hawkes$alpha*diag(est_hawkes$C%*%t(Gderivs[[2]])))-2*sum(diag(est_hawkes$C%*%t(Aderivs[[2]]))))/(p*T)

  ## Derivative with respect to gamma and alpha
  Sigma[q+1,(q+2):(q+1+p)] <- 2*diag(est_hawkes$C%*%t(Gderivs[[1]]))/(p*T)
  Sigma[(q+2):(q+1+p),q+1] <- Sigma[q+1,(q+2):(q+1+p)]

  ## Derivative with respect to gamma and C
  Sigma[q+1,(q+1+p+1):(q+1+p+p^2)] <- as.numeric(est_hawkes$C%*%Gammaderivs[[1]]+est_hawkes$C%*%t(Gammaderivs[[1]])+2*Matrix::Diagonal(p,est_hawkes$alpha)%*%Gderivs[[1]]-2*Aderivs[[1]])/(p*T)
  Sigma[(q+1+p+1):(q+1+p+p^2),q+1] <- Sigma[q+1,(q+1+p+1):(q+1+p+p^2)]

  ## Second derivative with respect to alpha
  Sigma[(q+2):(q+1+p),(q+2):(q+1+p)] <- 2*diag(V_diag)/(p*T)

  ## Compute derivative with respect to alpha and C
  Sigma[(q+2):(q+1+p),(q+1+p+1):(q+1+p+p^2)] <- 2*matrix(.Call("compute_alphaC_deriv",G),nrow=p)/(p*T)
  Sigma[(q+1+p+1):(q+1+p+p^2),(q+2):(q+1+p)] <- t(Sigma[(q+2):(q+1+p),(q+1+p+1):(q+1+p+p^2)])

  ## Compute Second Derivative with respect to C
  Sigma[(q+1+p+1):(q+1+p+p^2),(q+1+p+1):(q+1+p+p^2)] <- matrix(.Call("compute_d2C_deriv",Gamma),nrow=p^2)/(p*T)

  #### Compute Nodewise LASSO for the first columns of Sigma corresponding to beta and gamma
  ## Compute Nodewise LASSO using sigma from the paper
  Theta_tilde <- matrix(NA,ncol=1+q+p+p^2,nrow=q+1)
  sparsity <- 1+sum(est_hawkes$C!=0)/p+sum(rowSums(est_hawkes$C)^2)/p
  pen_weight <- 1/(p^(3/2)*log(p*T)^4*sparsity)
  for(j in 1:(q+1)) {
    Z <- as.numeric(tildeX%*%Sigma[,j])
    M <- as.matrix(tildeX%*%Sigma[,-j])
    m <- length(Z)
    node_lasso_Zsd <- sd(Z)*sqrt((m-1)/m)
    node_lasso_Msd <- apply(M,2,sd)*sqrt((m-1)/m)
    nvars <- dim(Sigma)[2]-1

    weight_scaling <- nvars/sum(1/node_lasso_Msd)

    node_wise_lasso <- glmnet::glmnet(t(t(M)/node_lasso_Msd),Z/node_lasso_Zsd,intercept=FALSE,standardize=FALSE,penalty.factor=weight_scaling/node_lasso_Msd,thresh=debias_thresh,maxit=debias_maxit)
    vec <- node_lasso_Zsd*coef(node_wise_lasso,s=pen_weight/(node_lasso_Zsd*m*weight_scaling),exact=debias_exact,x=t(t(M)/node_lasso_Msd),y=Z/node_lasso_Zsd,intercept=FALSE,standardize=FALSE,penalty.factor=weight_scaling/node_lasso_Msd,thresh=debias_thresh,maxit=debias_maxit)[-1]/node_lasso_Msd

    tau <- as.numeric((Sigma%*%Sigma)[j,j]-matrix((Sigma%*%Sigma)[j,-j],nrow=1)%*%vec)
    if(tau==0) {
      ## If tau=0 replace by one and print a warning, this is not ideal
      warning("The Sigma matrix in the debiasing has an almost zero column.")
      tau <- 1
    }

    Theta_tilde[j,-j] <- -vec/tau
    Theta_tilde[j,j] <- 1/tau
  }
  Theta <- Theta_tilde%*%Sigma

  ## Compute De-Biased Estimator
  theta_debiased <- c(est_hawkes$beta,est_hawkes$gamma)-Theta%*%matrix(grad,ncol=1)

  return(list(grad=grad,Sigma=Sigma,Theta=Theta,beta_debiased=theta_debiased[1:q],gamma_debiased=theta_debiased[q+1],tuning_parameter=pen_weight))
}

#' Computes a Complete Estimator with All Stages
#'
#' `NetHawkes` combines the functions `estimate_hawkes` and `debias_Hawkes` in a
#' two-step procedure: Estimate all parameters in the first stage, de-bias
#' `beta` and `gamma`, and fit, in a second stage, `C` and `alpha` keeping the
#' de-biased values of `beta` and `gamma` fixed.
#'
#' @inheritParams estimate_hawkes
#' @inheritParams debias_Hawkes
#' @param print.level Passed to [estimate_hawkes()]. If positive, in addition,
#'   information about in which stage the estimation is, is printed. The default
#'   is 0.
#' @param observation_matrix_network,observation_matrix_debiasing Similarly as
#'   in [estimate_hawkes()] these matrices are automatically computed when
#'   `NULL` is provided here (the default). Since the matrix is the same in
#'   repeated calls, it can lead to a speed up to compute the matrices once
#'   using [create_observation_matrix()]. The `observation_matrix_network` is of
#'   dimension `p` and the `observation_matrix_debiasing` is of dimension
#'   `1+q+p+p^2-1`. Here `p` denotes the number of vertices and `q` the
#'   dimension of the covariates.
#'
#' @return `NetHawkes` returns a list containing the elements:
#'   * `first_stage`: The estimator from the first stage as returned by
#'     [estimate_hawkes()]. It is computed using no initial values and
#'     estimating all model parameters.
#'   * `second_stage`: The estimator from the second stage as returned by
#'     [estimate_hawkes()]. It is computed fixing the values of `beta` and
#'     `gamma` to the de-biased estimators, and updating only `C` and `alpha`.
#'   * `debiasing`: The output of [debias_Hawkes()] from de-biasing the first
#'     stage estimator.
#'
#' @export
NetHawkes <- function(covariates,hawkes,omega,omega_alpha,lb,ub,C.ind.pen=NULL,print.level=0,max_iteration=100,tol=0.00001,link=exp,observation_matrix_network=NULL,observation_matrix_debiasing=NULL,cluster=NULL,debias_thresh=1e-14,debias_maxit=100000000,debias_exact=TRUE) {
  ## Perform first stage estimation
  if(print.level>0) {
    cat("Perform the first stage estimation.\n")
  }
  est_first_stage <- estimate_hawkes(covariates=covariates,hawkes=hawkes,omega=omega,omega_alpha=omega_alpha,lb=lb,ub=ub,C.ind.pen=C.ind.pen,fit_theta=TRUE,print.level=print.level,max_iteration=max_iteration,tol=tol,beta_init=NULL,gamma_init=NULL,alpha_init=NULL,link=link,observation_matrix=observation_matrix_network,cluster=cluster)

  ## Debiasing
  if(print.level>0) {
    cat("Debias the first stage estimator.\n")
  }
  debiased_est <- debias_Hawkes(covariates=covariates,hawkes=hawkes,est_hawkes=est_first_stage,link=link,observation_matrix=observation_matrix_debiasing,debias_thresh=debias_thresh,debias_maxit=debias_maxit,debias_exact=debias_exact)

  ## Compute Network estimate with debiased estimator
  if(print.level>0) {
    cat("Compute second stage estimator.\n")
  }
  est_second_stage <- estimate_hawkes(covariates=covariates,hawkes=hawkes,omega=omega,omega_alpha=omega_alpha,lb=lb,ub=ub,C.ind.pen=C.ind.pen,fit_theta=FALSE,print.level=print.level,max_iteration=max_iteration,tol=tol,beta_init=debiased_est$beta_debiased,gamma_init=debiased_est$gamma_debiased,alpha_init=est_first_stage$alpha,link=link,observation_matrix=observation_matrix_network,cluster=cluster)

  return(list(first_stage=est_first_stage,second_stage=est_second_stage,debiasing=debiased_est))
}



#' Compute a Complete Estimator Using Multiple Starting Points
#'
#' [NetHawkes_robust()] works similarly as [NetHawkes()] but it uses a
#' different optimization strategy. See Details for a precise description
#' of the differences.
#'
#' [NetHawkes_robust()] starts the optimization from `K` random starting
#' values for beta and gamma. If the user supplies some starting values in
#' `starting_beta` and `starting_gamma`, these will be extend by random values
#' selected uniformly between the supplied boundaries in `lb` and `ub`. Then,
#' the `nloptr` is used to optimize with respect to (beta,gamma), while for each
#' value of (beta,gamma), [estimate_hawkes()] will be run for fixed (beta,gamma)
#' to find the optimal value. Then, we take the minimum of all `K` results. These
#' runds of `nloptr` are computed with tolerance level given by 1000 times `tol`.
#' The optimal value is then refined by another run as before initialized with the
#' minimzer and tolerance level `tol`.
#'
#' All descrbibed changes affect the first stage estimator only. The later
#' stages remain the same as in [NetHawkes()]
#'
#' @inheritParams NetHawkes
#' @param K Number of starting values to use
#' @param starting_beta Matrix, each row of which corresponds to a starting
#'   value for beta to be used
#' @param starting_gamma Vector, each entry of which corresponds to a starting
#'   value for gamma to be used
#'
#' @return The retuned value is a list of the same structure as for [NetHawkes()]
#'   but with an additional element `nloptr` that contains the complete output
#'   of the refinement call from `nloptr`. This allows, e.g., to check for
#'   convergence of the optimization.
#'
#' @export
NetHawkes_robust <- function(covariates,hawkes,omega,omega_alpha,lb,ub,K,starting_beta=NULL,starting_gamma=NULL,C.ind.pen=NULL,print.level=0,max_iteration=100,tol=0.00001,link=exp,observation_matrix_network=NULL,observation_matrix_debiasing=NULL,cluster=NULL,debias_thresh=1e-14,debias_maxit=100000000,debias_exact=TRUE) {
  ## Read information
  q <- dim(covariates$cov[[1]])[2]

  ## Set information for optimisation
  args_init_opt <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=100*tol,print_level=0)
  args_refi_opt <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=tol,print_level=0)

  ## Create starting values
  starting_par <- matrix(NA,ncol=q+1,nrow=K)
  for(r in 1:q) {
    if(is.null(starting_beta)) {
      k0 <- 0
    } else {
      k0 <- dim(starting_beta)[1]
      starting_par[1:k0,r] <- starting_beta[,r]
    }
    if(k0<K) {
      starting_par[(k0+1):K,r] <- runif(K-k0,min=lb[r],max=ub[r])
    }
  }
  if(is.null(starting_gamma)) {
    k0 <- 0
  } else {
    k0 <- length(starting_gamma)
    starting_par[1:k0,q+1] <- starting_gamma
  }
  if(k0<K) {
    starting_par[(k0+1):K,q+1] <- runif(K-k0,min=lb[q+1],max=ub[q+1])
  }

  #### Robust first stage estimation
  ## Perform initial optimisation
  out <- list()
  obj_vals <- rep(NA,K)
  for(k in 1:K) {
    if(print.level>0) {
      cat("Initial estimation ",k," of ",K,".\n")
    }
    ## Check if provided starting value lies in the [lb,ub] range
    if(all(lb<=starting_par[k,]) & all(starting_par[k,]<=ub)) {
      intial_parameter_values <- starting_par[k,]
    } else {
      warning("In the initial step starting values outside the limits were produced, replace by average: This should not happen, contact the maintainer of this package")
      intial_parameter_values <- (lb+ub)/2
    }
    out[[k]] <- nloptr::nloptr(intial_parameter_values,estimate_hawkes_theta_container,opts=args_init_opt,ub=ub,lb=lb,covariates=covariates,hawkes=hawkes,omega=omega,omega_alpha=omega_alpha,C.ind.pen=C.ind.pen,print.level=print.level,max_iteration=max_iteration,tol=tol,alpha_init=NULL,link=link,observation_matrix=observation_matrix_network,cluster=cluster)
    obj_vals[k] <- out[[k]]$objective
  }

  ## Find minimum
  k0 <- min(which(obj_vals==min(obj_vals)))

  ## Run refining optimisation from optimal value
  if(print.level>0) {
    cat("Refinement step\n")
  }

  ## Check if computed value lies in the [lb,ub] region
  if(all(lb<=out[[k0]]$solution) & all(out[[k0]]$solution<=ub)) {
    ## It does
    intial_parameter_values_refinement <- out[[k0]]$solution
  } else {
    ## It does not, for some reason this can happen, replace the violating values
    intial_parameter_values_refinement <- out[[k0]]$solution
    if(!all(lb<=out[[k0]]$solution)) {
      lower_bound_violators <- which(lb>out[[k0]]$solution)
      intial_parameter_values_refinement[lower_bound_violators] <- lb[lower_bound_violators]+(ub[lower_bound_violators]-lb[lower_bound_violators])/100
    }
    if(!all(out[[k0]]$solution<=ub)) {
      upper_bound_violators <- which(ub<out[[k0]]$solution)
      intial_parameter_values_refinement[upper_bound_violators] <- ub[upper_bound_violators]-(ub[upper_bound_violators]-lb[upper_bound_violators])/100
    }

    print("In the refinement step starting values outside the limits were produced, replace by bound")

    print("Lower bound:")
    print(lb)
    print("Upper bound:")
    print(ub)
    print("Actual value")
    print(out[[k0]]$solution)
    print("Lower bound compare")
    print(lb<=out[[k0]]$solution)
    print("Upper bound compare")
    print(out[[k0]]$solution<=ub)
    print("Value to use")
    print(intial_parameter_values_refinement)
  }
  refined_out <- nloptr::nloptr(intial_parameter_values_refinement,estimate_hawkes_theta_container,opts=args_refi_opt,ub=ub,lb=lb,covariates=covariates,hawkes=hawkes,omega=omega,omega_alpha=omega_alpha,C.ind.pen=C.ind.pen,print.level=print.level,max_iteration=max_iteration,tol=tol,alpha_init=NULL,link=link,observation_matrix=observation_matrix_network,cluster=cluster)

  ## Run last estimate_hawkes to obtain estimates for alpha and C.
  if(print.level>0) {
    cat("Run estimate_hawkes on optimal parameter\n")
  }
  eh_out <- estimate_hawkes(covariates=covariates,hawkes=hawkes,omega=omega,omega_alpha=omega_alpha,lb=NULL,ub=NULL,C.ind.pen=C.ind.pen,fit_theta=FALSE,print.level=print.level,max_iteration=max_iteration,tol=tol,beta_init=refined_out$solution[1:q],gamma_init=refined_out$solution[q+1],alpha_init=NULL,link=link,observation_matrix=observation_matrix_network,cluster=cluster)


  #### Debiasing
  if(print.level>0) {
    cat("Debias the first stage estimator.\n")
  }
  debiased_est <- debias_Hawkes(covariates=covariates,hawkes=hawkes,est_hawkes=eh_out,link=link,observation_matrix=observation_matrix_debiasing,debias_thresh=debias_thresh,debias_maxit=debias_maxit,debias_exact=debias_exact)

  ## Compute Network estimate with debiased estimator
  if(print.level>0) {
    cat("Compute second stage estimator.\n")
  }
  est_second_stage <- estimate_hawkes(covariates=covariates,hawkes=hawkes,omega=omega,omega_alpha=omega_alpha,lb=NULL,ub=NULL,C.ind.pen=C.ind.pen,fit_theta=FALSE,print.level=print.level,max_iteration=max_iteration,tol=tol,beta_init=debiased_est$beta_debiased,gamma_init=debiased_est$gamma_debiased,alpha_init=eh_out$alpha,link=link,observation_matrix=observation_matrix_network,cluster=cluster)

  return(list(first_stage=eh_out,second_stage=est_second_stage,debiasing=debiased_est,nloptr=refined_out))
}

#' Compute Observation Matrix
#'
#' The observation matrix of dimension n required to handle the missing
#' intercept in the model is computed. For details please see our paper.
#'
#' The usage of this function is explained in the documentations for
#' [estimate_hawkes()] and [debias_Hawkes()]
#'
#' @param n Integer, the dimension of the required matrix.
#'
#' @return The observation matrix of the required dimension.
#'
#' @export
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
  M <- Matrix::sparseMatrix(i=i,j=j,x=x)

  ## Last row requires modification for odd n
  if(n %%2 ==1) {
    M[n,n] <- 1/sqrt(2)
    M[n,n+1] <- 0
    M[n,n+od] <- -1/sqrt(2)
  }

  return(Matrix::t(M))
}

#' Plot Intensity Functions and Hawkes Processes
#'
#' `plot_count_intensities` plots the provided Hawkes processes along with their
#' intensity functions. A plot window that provides enough space to plot a graph
#' for each vertex must be available before calling the function.
#'
#' @inheritParams simulate_hawkes
#' @param hawkes A Hawkes process in the form of a list with at least the
#'   elements `EL` and optionally `intensities` in the form as described in
#'   [simulate_hawkes()].
#' @param names An optional list of names that will be used as title for the
#'   individual plots. The vertex with number i in `hawkes$EL` will be labelled
#'   `names[i]`.
#'
#' @returns `plot_count_intensities` generates `p` plots (one for each vertex).
#'   The plots show the realized Hawkes processes and, if provided, the
#'   corresponding intensities. Before calling the function a plot window of the
#'   necessary size has to be created.
#' @export
plot_count_intensities <- function(hawkes,T,names=NULL) {
  ## Check if intensities is provided
  int_provided <- "intensities" %in% names(hawkes)
  names_provided <- !is.null(names)

  ## Read Information
  p <- max(hawkes$EL[,3])
  if(int_provided) {
    M <- max(hawkes$intensities)
  }

  ## Plot
  for(i in 1:p) {
    ## Extract events of process i
    ind <- which(hawkes$EL[,3]==i)

    ## Write plot title
    if(names_provided) {
      plot_title <- names[i]
    } else {
      plot_title <- sprintf("Vertex %d",i)
    }

    ## Create Plot
    Y <- length(ind)+1
    plot(0,0,type="n",xlim=c(0,T),ylim=c(0,Y),axes=FALSE,main=plot_title,xlab="Time",ylab="Events")
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

    ## Plot Intensity if required
    if(int_provided) {
      lines(hawkes$EL[,4],hawkes$intensities[i,]/M*Y,lty=2)
      label_positions <- c(0,round((1:4)/5*Y),Y)
      axis(4,at=label_positions,labels=round(label_positions/Y*M,digits=2))
    }
  }
}

#' Visualize Estimated Hawkes Causal Model in a Graph
#'
#' An estimate returned via [estimate_hawkes()] is visualized as a network.
#' Moreover, the network is returned. The use of this functions requires the
#' package igraph.
#'
#' `plot_interactions` generates a weighted network. The edge weights correspond
#' to the estimated values in `C` (potentially scaled by the value of
#' `edge.scaling`). The sizes of the vertices are such that there area is
#' proportional to `(0.1+sqrt(alpha))^2`, where `alpha` is the baseline activity
#' of the corresponding vertex (potentially multiplied with `vertex.scaling`).
#' The vertex names in the graph can be specified through `vertex.names`.
#'
#' @param estHawkes An estimated Hawkes process as returned from
#'   [estimate_hawkes()].
#' @param vertex.scaling,edge.scaling Non-negative numbers (1, by default) which
#'   scale edge weights and vertex sizes.
#' @param vertex.names Potential names of the vertices in the output graph, can
#'   be `NULL` (the default).
#' @param show.plot Logical value, if `TRUE` (the default) the network is
#'   plotted, otherwise the network is only returned but not plotted.
#' @param ... Additional arguments passed to the plotting routine of igraph if
#'   `show.plot=TRUE`.
#'
#' @returns `plot_interactions` returns an igraph object with the properties
#'   described under details. If `show.plot=TRUE`, the network is in addition
#'   plotted.
#' @export
plot_interactions <- function(estHawkes,vertex.scaling=1,edge.scaling=1,vertex.names=NULL,show.plot=TRUE,...) {
  ## Create Network
  G <- igraph::graph_from_adjacency_matrix(edge.scaling*estHawkes$C,mode="directed",weighted="weight")

  ## Add alpha as vertex attribute
  igraph::V(G)$alpha <- estHawkes$alpha

  ## Change Vertex names if applicable
  if(!is.null(vertex.names)) {
    igraph::V(G)$name <- vertex.names
  }

  ## Create Plot if asked
  if(show.plot) {
    ## Generate Node Sizes
    node_sizes <- vertex.scaling*(0.1+sqrt(estHawkes$alpha))

    ## Plot
    igraph::plot.igraph(G,vertex.size=node_sizes,edge.width=igraph::E(G)$weight,...)
  }

  return(G)
}

#' Compute the tuning parameter omega according to theory
#'
#' `compute_omega` computes a vector of tuning parameters according to the
#' theoretic results provided in our paper.
#'
#' This function computes the tuning parameter as suggested through Lemma 3.5 in
#' our paper. However, we do not use N_0 as in the paper but compute the
#' stochastic integral from the proof. `gamma_bar` is the value for gamma that
#' will be used for this. The remaining parameterst have the same meaning as in
#' Lemma 3.5.
#'
#' @inheritParams estimate_hawkes
#' @param p The number of vertices in the network
#' @param T The end of the observation period
#' @param alpha3,gamma_bar,mu Additional parameters that are required for the
#'   computation, cf. Details.
#'
#' @returns `compute_omega` returns a vector of length `p` that contains the
#'   penalty parameter for each vertex. It can be, e.g., provided to
#'   estimate_hawkes() as `omega`.
#'
#' @export
compute_omega <- function(hawkes,p,T,alpha3,gamma_bar,mu=log(2)) {
  ## Test computation
  B <- 0
  for(j in 1:p) {
    times <- sort(hawkes$EL[hawkes$EL[,3]==j,4])
    if(length(times)>0) {
      for(k in 1:length(times)) {
        B <- max(c(B,4*sum(exp(-gamma_bar*(times[k]-times[1:k])))/T))
      }
    }
  }

  ## Compute phi
  phi_mu <- exp(mu)-mu-1

  ## Compute Integral
  int <- .Call("compute_Vd_int",as.integer(p),hawkes$EL,as.double(gamma_bar))

  ## Compute Vd
  Vd <- 16*mu*int/((mu-phi_mu)*T^2)+B^2*(log(p)+log(p*T)+alpha3*log(T))/(mu-phi_mu)

  ## Compute dn
  dn <- 2*sqrt(Vd*(log(p)+log(p*T)+alpha3*log(T)))+B*(log(p)+log(p*T)+alpha3*log(T))/3

  return(dn)
}

#' Estimate the parameters from a multiple Hawkes process
#'
#' `MultiHawkes_robust` estimates the parameters of a Causal Hawkes process as
#' described in the paper when several realisations of this process are provided.
#'
#' Suppose that you have several observation from the same Hawkes causal model
#' as in our paper, that is, the parameters are the same for all realisations but
#' the covariates might be potentially different. `MultiHawkes_robust` estimates
#' this parameter by optimizing the following loss function:
#'
#' Sum over k ( Least squares loss for realisation k) + Penalty
#'
#' where the penalty equals
#'
#' Sum over i( omega[i]*sum(C[i,])) + omega_alpha*sum(alpha).
#'
#' The estimation method is the same as in [NetHawkes_robust()]. De-biasing is
#' not supported for multi observations.
#'
#' @inheritParams NetHawkes_robust
#' @param multi_covariates A list each entry of which contains a `covariate` as
#'   described in, e.g., [NetHawkes_robust()]. Each entry of the list corresponds to
#'   another realisation of the Hawkes process.
#' @param multi_hawkes A list each entry of which contains a Hawkes process in
#'   the same format as for, e.g, [NetHawkes_robust()]. Each entry of
#'   `multi_hawkes` corresponds to one realisation of the Hawkes process.
#'
#' @returns `MultiHawkes_robust` returns a list containing the estimates for the
#'   parameters `C`, `alpha`, `beta`, and `gamma`.
#'
#' @export
MultiHawkes_robust <- function(multi_covariates,multi_hawkes,omega,omega_alpha,lb,ub,K,starting_beta=NULL,starting_gamma=NULL,C.ind.pen=NULL,print.level=0,max_iteration=100,tol=0.00001,alpha_init=NULL,link=exp,observation_matrix=NULL,cluster=NULL) {
  ## Read information
  q <- dim(multi_covariates[[1]]$cov[[1]])[2]

  ## Set information for optimisation
  args_init_opt <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=100*tol,print_level=0)
  args_refi_opt <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=tol,print_level=0)

  ## Create starting values
  starting_par <- matrix(NA,ncol=q+1,nrow=K)
  for(r in 1:q) {
    if(is.null(starting_beta)) {
      k0 <- 0
    } else {
      k0 <- dim(starting_beta)[1]
      starting_par[1:k0,r] <- starting_beta[,r]
    }
    if(k0<K) {
      starting_par[(k0+1):K,r] <- runif(K-k0,min=lb[r],max=ub[r])
    }
  }
  if(is.null(starting_gamma)) {
    k0 <- 0
  } else {
    k0 <- length(starting_gamma)
    starting_par[1:k0,q+1] <- starting_gamma
  }
  if(k0<K) {
    starting_par[(k0+1):K,q+1] <- runif(K-k0,min=lb[q+1],max=ub[q+1])
  }

  ## Perform initial optimisation
  out <- list()
  obj_vals <- rep(NA,K)
  for(k in 1:K) {
    if(print.level>0) {
      cat("Initial estimation ",k," of ",K,".\n")
    }
    out[[k]] <- nloptr::nloptr(starting_par[k,],estimate_theta_multi_hawkes,opts=args_init_opt,ub=ub,lb=lb,multi_covariates=multi_covariates,multi_hawkes=multi_hawkes,omega=omega,omega_alpha=omega_alpha,C.ind.pen=C.ind.pen,print.level=print.level,max_iteration=max_iteration,tol=tol,alpha_init=alpha_init,link=link,observation_matrix=observation_matrix,cluster=cluster,return_objective=TRUE)
    obj_vals[k] <- out[[k]]$objective
  }
  if(print.level>0) {
    cat("Found values of the objective:\n")
    print(obj_vals)
    cat("Corresponding parameters are:\n")
    for(k in 1:K) {
      print(out[[k]]$solution)
    }
  }

  ## Find minimum
  k0 <- min(which(obj_vals==min(obj_vals)))

  ## Run refining optimisation from optimal value
  if(print.level>0) {
    cat("Refinement step\n")
  }
  refined_out <- nloptr::nloptr(out[[k0]]$solution,estimate_theta_multi_hawkes,opts=args_refi_opt,ub=ub,lb=lb,multi_covariates=multi_covariates,multi_hawkes=multi_hawkes,omega=omega,omega_alpha=omega_alpha,C.ind.pen=C.ind.pen,print.level=print.level,max_iteration=max_iteration,tol=tol,alpha_init=alpha_init,link=link,observation_matrix=observation_matrix,cluster=cluster,return_objective=TRUE)

  ## Run last estimate_hawkes to obtain estimates for alpha and C.
  if(print.level>0) {
    cat("Run estimate_hawkes on optimal parameter\n")
  }
  eh_out <- estimate_theta_multi_hawkes(theta=refined_out$solution,multi_covariates=multi_covariates,multi_hawkes=multi_hawkes,omega=omega,omega_alpha=omega_alpha,C.ind.pen=C.ind.pen,print.level=print.level,max_iteration=max_iteration,tol=tol,alpha_init=alpha_init,link=link,observation_matrix=observation_matrix,cluster=cluster)

  return(eh_out)
}

#' Use cross-validation to find the tuning parameter for penalizing `C`
#'
#' `cvMultiHawkes` performs cross validation to find the tuning parameter to be
#' used for the estimation in [MultiHawkes_robust()].
#'
#' In order to perform cross-validation, the data is split in time in a training
#' data-set and a testing data-set. The time length of the training data-set is
#' `tf`*`T`, where `T` is the length of the observation period. The fitting
#' method is [MultiHawkes_robust()]. In order to evaluate the performance of a
#' given tuning parameter, the least squares loss is computed on the remaining
#' data, the testing data-set. The method to search for the optimal tuning
#' parameter is Golden-Section search with lower bound 0 and upper bound
#' provided in `omega_start`. While [MultiHawkes_robust()] is called on the
#' entire process, the least squares losses are computed for each vertex
#' separately. Thus, also the optimisation of the tuning parameter is performed
#' separately.
#'
#' The Golden-Section search is executed `M` times (regardless of convergence of
#' the algorithm). Finally, the midpoint of the resulting interval is evaluated.
#' The return value is the choice of the tuning parameter for each vertex that
#' yields the lowest individual loss.
#'
#' This routine is also suitable for cross-validation of a single realisation
#' Hawkes process whenn `multi_hawkes` and `multi_covariates` contain only one
#' element.
#'
#' @inheritParams MultiHawkes_robust
#' @param multi_covariates A list of the same format as in
#'   [MultiHawkes_robust()]. However, it is required that the `times` vector in
#'   each entry is the same.
#' @param omega_start A vector (its length must equal the number of vertices)
#'   that contains the starting values of the search. It should contain values
#'   that are too large (see also below for details)
#' @param nos The number of starting points to be used in each call of
#'   [MultiHawkes_robust()], it defaults to 5.
#' @param tf The fraction of the data (in time) that is used for training,
#'   defaults to 0.8.
#' @param M The number of steps to be used in the Golden-Section search,
#'   defaults to 5
#'
#' @returns A list containing the following three elements:
#'   *  `omega`: The optimal choice of the tuning parameter among all evaluated
#'               choices. This can be passed to, e.g, to [MultiHawkes_robust()]
#'               or [NetHawkes_robust()].
#'   * `computed_omega`: A matrix with 2M+1 rows. Each row contains a set of
#'                       tuning parameters that was evaluated.
#'   * `LSvals`: A matrix with 2M+1 rows. Each row contains the individual
#'              (per vertex) least squares losses of the corresponding tuning
#'              parameters provided in entry `computed_omega`.
#'
#' @export
cvMultiHawkes <- function(multi_Hawkes,multi_covariates,omega_start,lb,ub,nos=5,tf=0.8,M=5,starting_beta=NULL,starting_gamma=NULL,C.ind.pen=NULL,print.level=0,max_iteration=100,tol=0.00001,alpha_init=NULL,link=exp,observation_matrix=NULL,cluster=NULL) {
  ## Read information
  K <- length(multi_Hawkes)
  T <- max(multi_covariates[[1]]$times)
  p <- dim(multi_covariates[[1]]$cov[[1]])[1]
  phi <- (1+sqrt(5))/2

  ## Compute training time
  training_time <- max(multi_covariates[[1]]$times[multi_covariates[[1]]$times<tf*T])

  ## Split data in training and test data
  training_hawkes <- list()
  training_covariates <- list()
  test_hawkes <- list()
  test_covariates <- list()

  for(k in 1:K) {
    ## Find training and test indices in edge list
    trai_ind <- which(multi_Hawkes[[k]]$EL[,4]<=training_time)
    test_ind <- which(multi_Hawkes[[k]]$EL[,4]> training_time)
    trai_EL <- multi_Hawkes[[k]]$EL[trai_ind,]
    test_EL <- multi_Hawkes[[k]]$EL[test_ind,]
    test_EL[,4] <- test_EL[,4]-training_time

    ## Save Hawkes training and test Data
    training_hawkes[[k]] <- list(EL=trai_EL[order(trai_EL[,4]),])
    test_hawkes[[k]] <- list(EL=test_EL[order(test_EL[,4]),])

    ## Get training and test time indices in covariates
    trai_cov_ind <- which(multi_covariates[[k]]$times<=training_time)
    test_cov_ind <- which(multi_covariates[[k]]$times >training_time)

    ## Write new time vectors
    training_covariates[[k]] <- list(times=multi_covariates[[k]]$times[trai_cov_ind],cov=list())
    test_covariates[[k]] <- list(times=multi_covariates[[k]]$times[test_cov_ind]-training_time,cov=list())

    ## Copy covariates to respective data-set
    last_training_index <- max(trai_cov_ind)
    for(i in 1:last_training_index) {
      ## Add to training set
      training_covariates[[k]]$cov[[i]] <- multi_covariates[[k]]$cov[[i]]
    }
    for(i in (last_training_index+1):length(multi_covariates[[k]]$cov)) {
      ## Add to test set
      test_covariates[[k]]$cov[[i-last_training_index]] <- multi_covariates[[k]]$cov[[i]]
    }
  }

  ## Prepare result matrices
  computed_omega <- matrix(NA,nrow=2*M+1,ncol=p)
  LSvals <- matrix(0,nrow=2*M+1,ncol=p)
  omega_ub <- omega_start
  omega_lb <- rep(0,p)

  ## Perform Golden-Section search
  for(m in 1:M) {
    cat("Run ",m," of ",M,".\n")
    ## Compute interior points
    omega1 <- omega_ub-(omega_ub-omega_lb)/phi
    omega2 <- omega_lb+(omega_ub-omega_lb)/phi

    ## Save omegas
    computed_omega[2*m-1,] <- omega1
    computed_omega[2*m  ,] <- omega2

    ## Compute estimates
    est1 <- MultiHawkes_robust(training_covariates,training_hawkes,omega1,omega_alpha=0,lb=lb,ub=ub,K=nos,starting_beta=NULL,starting_gamma=NULL,C.ind.pen=NULL,print.level=print.level,max_iteration=max_iteration,tol=tol,alpha_init=alpha_init,link=link,observation_matrix=observation_matrix,cluster=cluster)
    est2 <- MultiHawkes_robust(training_covariates,training_hawkes,omega2,omega_alpha=0,lb=lb,ub=ub,K=nos,starting_beta=NULL,starting_gamma=NULL,C.ind.pen=NULL,print.level=print.level,max_iteration=max_iteration,tol=tol,alpha_init=alpha_init,link=link,observation_matrix=observation_matrix,cluster=cluster)

    ## Compute corresponding least squares
    for(k in 1:K) {
      LSvals[2*m-1,] <- LSvals[2*m-1,]+compute_individual_lest_squares_theta(c(est1$beta,est1$gamma),test_covariates[[k]],est1$C,est1$alpha,test_hawkes[[k]],link=link)
      LSvals[2*m  ,] <- LSvals[2*m  ,]+compute_individual_lest_squares_theta(c(est2$beta,est2$gamma),test_covariates[[k]],est2$C,est2$alpha,test_hawkes[[k]],link=link)
    }

    ## Compute new bounds
    if(sum(est1$C)==0) {
      ## Lower choice of omega yields empty network -> Set omega1 as upper bound
      omega_ub <- omega1
    } else {
      ## Existent network in estimation, use Golden-Section update
      upper_better <- which(LSvals[2*m-1,]> LSvals[2*m,])
      lower_better <- which(LSvals[2*m-1,]<=LSvals[2*m,])
      omega_lb[upper_better] <- omega1[upper_better]
      omega_ub[lower_better] <- omega2[lower_better]
    }
  }

  ## Compute Estimate for the mid point of the resulting interval
  omega_mid <- (omega_lb+omega_ub)/2
  computed_omega[2*M+1,] <- omega_mid
  est <- MultiHawkes_robust(training_covariates,training_hawkes,omega_mid,omega_alpha=0,lb=lb,ub=ub,K=nos,starting_beta=NULL,starting_gamma=NULL,C.ind.pen=NULL,print.level=print.level,max_iteration=max_iteration,tol=tol,alpha_init=alpha_init,link=link,observation_matrix=observation_matrix,cluster=cluster)
  for(k in 1:K) {
    LSvals[2*M+1,] <- LSvals[2*M+1,]+compute_individual_lest_squares_theta(c(est$beta,est$gamma),test_covariates[[k]],est$C,est$alpha,test_hawkes[[k]],link=link)
  }

  ## Find minimum
  return_omega <- rep(NA,p)
  for(i in 1:p) {
    min_ind <- which(LSvals[,i]==min(LSvals[,i]))
    return_omega[i] <- min(computed_omega[min_ind,i])
  }

  return(list(omega=return_omega,computed_omega=computed_omega,LSvals=LSvals))
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
  obj <- compute_lest_squares_theta(par=theta,covariates=covariates,C=opt_theta$C,alpha=opt_theta$alpha,hawkes=hawkes,link=link)/T+2*sum(omega*opt_theta$C)

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
    return(as.numeric((t(alpha)%*%V%*%alpha+sum(diag(C%*%Gamma))+2*sum(alpha*diag(C%*%t(G)))-2*sum(alpha*v)-2*sum(diag(C%*%t(A))))/(p*T)))
  } else {
    return(list(C=C,alpha=alpha,beta=beta,gamma=gamma))
  }

}
