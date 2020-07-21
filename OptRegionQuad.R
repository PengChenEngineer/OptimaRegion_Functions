OptRegionQuad <- function(design, y, constr_lb, constr_ub,
                          triangularRegion = FALSE, vertex1 = NULL, vertex2 = NULL,
                          alpha = 0.05, B = 200, maximization = TRUE) {
  X <- design
  nosim <- B
  LB <- constr_lb
  UB <- constr_ub
  # Check this is for k=2 factors only
  k <- dim(X)[2]
  if ((k > 2) | (k < 2)) stop("Error. Number of factors must equal to 2")
  
  # If experimental region was specified as triangular, compute the parameteres defining the 3 lines that approximate its shape.
  if (triangularRegion) {
    x11 <- vertex1[1] # use user defined vertices; third vertex is (0,0)
    x21 <- vertex1[2]
    x12 <- vertex2[1]
    x22 <- vertex2[2]
    m1 <- x21 / x11
    m2 <- x22 / x12
    m3 <- (x21 - x22) / (x11 - x12)
    bintercept <- x22 - m3 * x12
  } else {
    m1 <- 0
    m2 <- 0
    m3 <- 0
    bintercept <- 0
  }
  
  # Find thetaHat and V, the coefficients and variance-covariance matrix for a second order (quadratic) polynomial model
  data <- data.frame(cbind(X, y))
  colnames(data) <- c("x1", "x2", "y")
  # Full quadratic polynomial assumed
  model <- lm(y ~ rsm::FO(x1, x2) + rsm::TWI(x1, x2) + rsm::PQ(x1, x2), data = data)
  names(model$coefficients) <- NULL
  thetaHat <- model$coefficients
  V <- vcov(model)
  p <- length(thetaHat)
  numberN <- length(y)
  # Nonparametric bootstrapping of the residuals used to simulate Theta
  fit <- fitted(model)
  names(fit) <- NULL
  # res<-resid(model)
  # names(res)<-NULL
  res <- resid(model) - mean(resid(model)) # correct for bias
  # Use studentized residuals to correct for bias
  # res<-rstudent(model)-mean(rstudent(model))
  Theta <- matrix(nrow = nosim, ncol = p)
  ThetaStd <- matrix(nrow = nosim, ncol = p)
  print(c("Bootstrapping..."))
  for (i in 1:nosim) {
    indices <- sample(1:numberN, replace = TRUE)
    ystar <- fit + res[indices]
    datastar <- data.frame(cbind(X, ystar))
    colnames(datastar) <- c("x1", "x2", "ystar")
    # Full quadratic polynomial assumed
    modelstar <- lm(ystar ~ rsm::FO(x1, x2) + rsm::TWI(x1, x2) + rsm::PQ(x1, x2), data = datastar)
    names(modelstar$coefficients) <- NULL
    Theta[i, ] <- modelstar$coefficients
    Vi <- vcov(modelstar)
    # find square root matrix and standardize Theta vectors following Yeh and Singh, J.Royal Stat. Soc. 59(3), pp. 639-652, 1997
    e <- eigen(Vi)
    Ssqrtinv <- solve(e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)) # use solve, not ^(1)
    ThetaStd[i, ] <- Ssqrtinv %*% (Theta[i, ] - thetaHat) * sqrt(numberN)
  }
  print(c("Computing parameter vectors inside their CR and finding maxima..."))
  d <- vector(length = nosim)
  ThetavectorsStdMat <- as.matrix(ThetaStd)
  # Calculate data depth of each parameter vector
  d <- DepthProc::depthTukey(ThetavectorsStdMat, ThetavectorsStdMat, ndir = 3000)
  order.d <- order(d) # order only based on Tukey's depth
  # find number of points should be in the alpha percent confidence region
  ind.alpha <- alpha * nosim + 1
  # provide the indices of points in the confidence region
  indices <- order.d[(ind.alpha:nosim)]
  # Compute vector with Theta values inside the CR
  Theta_In <- Theta[indices, ]
  l <- dim(Theta_In)
  # First create initial points for optimizer
  delta1 <- 0.1 * (UB[1] - LB[1])
  delta2 <- 0.1 * (UB[2] - LB[2])
  X0 <- matrix(nrow = 5, ncol = k)
  if (triangularRegion) {
    X0[1, ] <- c(LB[1] + delta1, LB[2] + delta2)
    X0[2, ] <- vertex1 - c(0, delta2)
    X0[3, ] <- vertex2 - c(delta1, 0)
    X0[4, ] <- c(mean(X0[1:3, 1]), mean(X0[1:3, 2]))
    # make sure points are inside LB and UB
    X0[, 2] <- apply(cbind(X0[, 2], UB[2]), 1, min)
    X0[, 2] <- apply(cbind(X0[, 2], LB[2]), 1, max)
    X0[, 1] <- apply(cbind(X0[, 1], UB[1]), 1, min)
    X0[, 1] <- apply(cbind(X0[, 1], LB[1]), 1, max)
    noinitial <- 4
    numberFeasible <- 4
  } else {
    X0[1, ] <- c(LB[1] + delta1, LB[2] + delta2)
    X0[2, ] <- c(LB[1] + delta1, UB[2] - delta2)
    X0[3, ] <- c(UB[1] - delta1, LB[2] + delta2)
    X0[4, ] <- c(UB[1] - delta1, UB[2] - delta2)
    X0[5, ] <- c(mean(X0[1:4, 1]), mean(X0[1:4, 2]))
    noinitial <- 5
    numberFeasible <- 5
  }
  xin <- matrix(nrow = l[1], ncol = k)
  # Main optimization loop
  for (m in 1:l[1]) {
    # Pick a Theta inside the alpha% data depth region
    betaCoef <- Theta_In[m, ]
    best <- 1e50
    for (j in 1:numberFeasible) {
      # Setup optimization problem depending on whether the experimental region has been defined as triangular or not
      if (triangularRegion) {
        out <- nloptr::nloptr(X0[j, ], eval_f = computefQuad, eval_grad_f = computeg, eval_jac_g_ineq = compute_Jacobian_Const, lb = LB, ub = UB, eval_g_eq = NULL, eval_g_ineq = constraintsQuad, opts = list("algorithm" = "NLOPT_LD_MMA", print_level = 0, xtol_rel = 1e-03), betaCoef = betaCoef, maximization = maximization, m1 = m1, m2 = m2, m3 = m3, bintercept = bintercept)
      } else {
        out <- nloptr::nloptr(X0[j, ], eval_f = computefQuad, eval_grad_f = computeg, lb = LB, ub = UB, eval_g_eq = NULL, eval_g_ineq = NULL, opts = list("algorithm" = "NLOPT_LD_MMA", print_level = 0, xtol_rel = 1e-03), betaCoef = betaCoef, maximization = maximization, m1 = m1, m2 = m2, m3 = m3, bintercept = bintercept)
      }
      if ((out$objective < best) & (out$status > 0)) {
        best <- out$objective
        bestSol <- out$sol
        bestStatus <- out$status
      }
      # print(c('in=',X0[j,],'out=',bestSol))
    } # endfor j
    # save best solution found among all tries for simulated parameter set m
    xin[m, ] <- bestSol
    # print(c(m, best, xin[m, ], bestStatus))
  } # endfor m
  
  # return
  optima <- data.frame(xin)
  opt_bag <- apply(xin, 2, mean)
  opt_bag <- data.frame(X1 = opt_bag[1], X2 = opt_bag[2])
  design <- data.frame(X1 = design[, 1], X2 = design[, 2])
  structure(
    list(
      optima = optima, opt_bag = opt_bag,
      design = design, y = y,
      constr_lb = constr_lb, constr_ub = constr_ub
    ),
    class = "crquad"
  )
} 

# -------------------------------------------------------------------------

plot.crquad <- function(res, xlab = "x1", ylab = "x2") {
  xin <- res$optima
  centroid <- res$opt_bag
  X <- res$design
  y <- res$y
  LB <- res$constr_lb
  UB <- res$constr_ub
  tpsfit <- fields::Tps(X, y, lambda = 0.04)
  surface <- fields::predictSurface(tpsfit)
  image(surface, lwd = 2, col = heat.colors(0), xlim = c(LB[1], UB[1]), ylim = c(LB[2], UB[2]),
        xlab = xlab, ylab = ylab)
  contour(surface, add = T, drawlabels = T, lwd = 2,  xlim = c(LB[1], UB[1]), ylim = c(LB[2], UB[2]))
  hpts_original <- chull(xin)
  hpts_closed <- c(hpts_original, hpts_original[1])
  lines(xin[hpts_closed, ], col = "blue")
  polygon(xin[hpts_closed, 1], xin[hpts_closed, 2], col = "grey")
  points(centroid[1], centroid[2], col = "red", pch = 19)
}

# -------------------------------------------------------------------------

constraintsQuad <- function(x, betaCoef, maximization, m1, m2, m3, bintercept) {
  # Computes the constraints limiting the triangular-like experimental region (approximated with a triangle, so 3 constraints)
  z <- vector(length = 3)
  z[1] <- x[2] - m1 * x[1]
  z[2] <- m2 * x[1] - x[2]
  z[3] <- x[2] - bintercept - m3 * x[1]
  return(z)
}

# -------------------------------------------------------------------------

plotConvexHull <- function(xin, LB, UB, xlab = "X", ylab = "Y") {
  # Plots the convex hull of the points in vector xin
  plot(0, 0, col = "white", xlim = c(LB[1], UB[1]), ylim = c(LB[2], UB[2]), xlab = xlab, ylab = ylab)
  hpts_original <- chull(xin)
  hpts_closed <- c(hpts_original, hpts_original[1])
  lines(xin[hpts_closed, ], col = "blue")
  polygon(xin[hpts_closed, 1], xin[hpts_closed, 2], col = "grey")
  return(hpts_original)
}

# -------------------------------------------------------------------------

computefQuad <- function(x, betaCoef, maximization, m1, m2, m3, bintercept) {
  # Returns -f(x)--a scalar (second order quadratic polynomial in k variables, where k=2 or 3)
  k <- length(x)
  # assume model has an intercept; CONFORMS TO MATLAB X2FX order of coefficients
  if (k == 2) {
    f <- betaCoef[1] + x[1] * betaCoef[2] + x[2] * betaCoef[3] + x[1] * x[2] * betaCoef[4] + x[1]^2 * betaCoef[5] + x[2]^2 * betaCoef[6]
  } else {
    # k=3
    f <- betaCoef[1] + x[1] * betaCoef[2] + x[2] * betaCoef[3] + x[3] * betaCoef[4] + x[1] * x[2] * betaCoef[5] + x[1] * x[3] * betaCoef[6] + x[2] * x[3] * betaCoef[7] + x[1]^2 * betaCoef[8] + x[2]^2 + betaCoef[9] + x[3]^2 * betaCoef[10]
  }
  if (maximization) {
    return(-f)
  } else {
    return(f)
  }
}

# -------------------------------------------------------------------------

computeg <- function(x, betaCoef, maximization, m1, m2, m3, bintercept) {
  # Returns gradient of f(x) evaluated at x (a vector)
  k <- length(x)
  g <- vector(length = k)
  # assume model has an intercept; CONFORMS TO MATLAB X2FX
  if (k == 2) {
    g[1] <- betaCoef[2] + 2 * betaCoef[5] * x[1] + betaCoef[4] * x[2]
    g[2] <- betaCoef[3] + 2 * betaCoef[6] * x[2] + betaCoef[4] * x[1]
  } else {
    # k=3
    g[1] <- betaCoef[2] + 2 * betaCoef[8] * x[1] + betaCoef[5] * x[2] + betaCoef[6] * x[3]
    g[2] <- betaCoef[3] + 2 * betaCoef[9] * x[2] + betaCoef[5] * x[1] + betaCoef[7] * x[3]
    g[3] <- betaCoef[4] + 2 * betaCoef[10] * x[3] + betaCoef[6] * x[1] + betaCoef[7] * x[2]
  }
  if (maximization) {
    return(-g)
  } else {
    return(g)
  }
}

# -------------------------------------------------------------------------

compute_Jacobian_Const <- function(x, betaCoef, maximization, m1, m2, m3, bintercept) {
  return(rbind(
    c(-m1, 1),
    c(m2, -1),
    c(-m3, 1)
  ))
}

# -------------------------------------------------------------------------

computeH <- function(x, betaCoef, maximization) {
  # Evaluates the Hessian at x
  k <- length(x)
  H <- matrix(nrow = k, ncol = k)
  if (k == 2) {
    H[1, 1] <- 2 * betaCoef[5]
    H[2, 2] <- 2 * betaCoef[6]
    H[1, 2] <- betaCoef[4]
    H[2, 1] <- H[1, 2]
  } else {
    # k=3
    H[1, 1] <- 2 * betaCoef[8]
    H[2, 2] <- 2 * betaCoef[9]
    H[3, 3] <- 2 * betaCoef[10]
    H[1, 2] <- betaCoef[5]
    H[1, 3] <- betaCoef[6]
    H[2, 3] <- betaCoef[7]
    H[2, 1] <- H[1, 2]
    H[3, 1] <- H[1, 3]
    H[3, 2] <- H[2, 3]
  }
  if (maximization) {
    return(-H)
  }
  else {
    return(H)
  }
}