GloptiPolyRegion <- function(design, y, constr_lb, constr_ub,
                             alpha = 0.05, B = 200, degree,  
                             maximization = TRUE, verbose = TRUE) {
  X <- design
  lb <- constr_lb
  ub <- constr_ub
  
  X <- data.frame(X)
  y <- data.frame(y)
  # Check polynomial order 
  if (degree < 2 || degree > 3) {
    stop("This function accepts only quadratic or cubic polynomials!")
  }
  # Check number of variables 
  if (ncol(X) < 2 || ncol(X) > 5) {
    stop("This function accepts only 2 - 5 variables!")
  }
  # Simplify function arguments 
  scale <- TRUE # scale X to [-1, 1]
  # Original fit 
  if (verbose) print("Fit original model ...")
  if (scale) X <- encode_orthogonal(X, lb, ub) # scale X to [-1, 1]
  model_original <- fit_polym(X, y, degree) # fit polynomial model
  theta_hat <- model_original$coefficients
  cov_theta_hat <- vcov(model_original)
  p <- length(theta_hat)
  n <- nrow(y)
  # Bootstrap 
  if (verbose) print("Bootstrap ...")
  theta_hat_star <- matrix(NA, nrow = B, ncol = p)
  theta_hat_star_studentized <- matrix(NA, nrow = B, ncol = p)
  boot_optima <- matrix(NA, nrow = B, ncol = ncol(X))
  fit <- fitted(model_original)
  names(fit) <- NULL
  res <- resid(model_original) - mean(resid(model_original))
  # accumulate B "optimizable" bootstrap instances
  counter <- 0
  while (counter < B) {
    model_star <- fit_polym(X, fit + res[sample(1:n, replace = TRUE)], degree)
    opt <- try(
      # optimize a polynomial model fitted through the polym formula using
      # the GloptiPolyR solver
      Gloptipolym(
        coefficients = model_star$coefficients,
        k = ncol(X), degree = degree,
        lb = rep(-1, ncol(X)), ub = rep(1, ncol(X)),
        maximization = maximization
      )
    )
    if (!inherits(opt, "try-error")) {
      counter <- counter + 1
      if (verbose) print(paste(counter, "th bootstrapping ..."))
      theta_hat_star[counter, ] <- model_star$coefficients
      e <- eigen(vcov(model_star))
      S_sqrt_inv <- solve(e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors))
      theta_hat_star_studentized[counter, ] <- sqrt(n) * S_sqrt_inv %*%
        (theta_hat_star[counter, ] - theta_hat)
      boot_optima[counter, ] <- opt$solution
    }
  }
  # Trim 
  if (verbose) print("Trimming ...")
  d <- vector(length = B)
  d <- DepthProc::depthTukey(theta_hat_star_studentized,
                             theta_hat_star_studentized,
                             ndir = 3000
  )
  order_d <- order(d)
  ind_alpha <- alpha * B + 1
  indices <- order_d[ind_alpha:B]
  # Optimization 
  # extract already optimized results that remains after trimming
  if (verbose) print("Optimizing over bootstrapped models that remains...")
  boot_optima <- boot_optima[indices, ]
  if (scale) boot_optima <- decode_orthogonal(boot_optima, lb, ub)
  
  # return
  optima <- data.frame(boot_optima)
  opt_bag <- apply(boot_optima, 2, mean)
  opt_bag <- data.frame(
    X1 = opt_bag[1], X2 = opt_bag[2], X3 = opt_bag[3], X4 = opt_bag[4], X5 = opt_bag[5]
  )
  structure(
    list(
      optima = optima, opt_bag = opt_bag,
      constr_lb = constr_lb, constr_ub = constr_ub
    ),
    class = "crpoly"
  )
}

# -------------------------------------------------------------------------

plot.crpoly <- function(res, axes_labels) {
  boot_optima <- res$optima
  bagged_optimum <- res$opt_bag
  lb <- res$constr_lb
  ub <- res$constr_ub
  draw_2D_CRs(
    boot_optima, bagged_optimum, lb, ub,
    axes_labels = axes_labels
  )
}

# -------------------------------------------------------------------------

draw_2D_CR <- function(boot_optima, boost_optimum,
                       xlab, ylab, xlim, ylim) {
  # get the indices of the points that are on the boundary of the convex hull
  id_cvx_hull <- chull(boot_optima)
  id_cvx_hull <- c(id_cvx_hull, id_cvx_hull[1]) # add 1st point to get a loop
  # cluster the bootstrap optima
  cluster_boot_optima <- try( # in case of perfect fit
    mclust::densityMclust(boot_optima, verbose = FALSE),
    silent = TRUE
  )
  if (!inherits(cluster_boot_optima, "try-error")) {
    # plot contours of estimated density of the bootstrap optima
    plot(cluster_boot_optima,
         what = "density", type = "hdr", col = "gray",
         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim
    )
    # plot the bootstrap optima points
    points(boot_optima,
           col = "black", cex = 0.5, pch = 16,
           xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim
    )
  } else {
    plot(boot_optima,
         col = "black", cex = 0.5, pch = 16,
         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim
    )
  }
  # plot the boundary of the convex hull of the bootstrap optima
  lines(boot_optima[id_cvx_hull, ], col = "black")
  # plot the boosted optimum
  points(boost_optimum[1], boost_optimum[2], col = "red", cex = 1, pch = 16)
}

# -------------------------------------------------------------------------

draw_2D_CRs <- function(boot_optima, boost_optimum, lb, ub,
                        for_dev = TRUE, axes_labels) {
  if (for_dev) dev.new()
  k <- ncol(boot_optima)
  par(mfrow = c(k - 1, k - 1))
  for (i in 1:(k - 1)) { # each row of the sub-figures
    for (j in 2:k) { # sub-figure in each row
      if (i < j) {
        if (is.null(axes_labels)) {
          xlab <- paste("x", j)
          ylab <- paste("x", i)
        } else {
          if(!is.null(axes_labels) && length(axes_labels) != k) {
            stop("Incorrect number of labels!")
          }
          xlab <- axes_labels[j]
          ylab <- axes_labels[i]
        }
        draw_2D_CR(
          boot_optima = boot_optima[, c(j, i)],
          boost_optimum = boost_optimum[c(j, i)],
          xlab = xlab, ylab = ylab,
          xlim = c(lb[j], ub[j]), ylim = c(lb[i], ub[i])
        )
      } else {
        plot.new()
      }
    }
  }
}

# -------------------------------------------------------------------------

encode_orthogonal <- function(Xi, lb, ub) {
  Xi <- t(Xi) # put data in columns for vectorized implementation
  M <- (lb + ub) / 2
  dim(M) <- c(nrow(Xi), 1)
  res <- sweep(Xi, 1, M, "-")
  R <- ub - lb
  S <- diag(R) / 2
  res <- solve(S, res)
  # return
  t(res)
}

# -------------------------------------------------------------------------

decode_orthogonal <- function(X, lb, ub) {
  X <- t(X)
  M <- (lb + ub) / 2
  dim(M) <- c(nrow(X), 1)
  R <- ub - lb
  S <- diag(R) / 2
  res <- S %*% X
  res <- sweep(res, 1, M, "+")
  # return
  t(res)
}

# -------------------------------------------------------------------------

fit_polym <- function(X, y, degree) {
  data <- data.frame(X, y)
  if (ncol(X) == 2) {
    colnames(data) <- c("x1", "x2", "y")
    model <- lm(y ~ polym(x1, x2, degree = degree, raw = TRUE), data = data)
  } else if (ncol(X) == 3) {
    colnames(data) <- c("x1", "x2", "x3", "y")
    model <- lm(y ~ polym(x1, x2, x3, degree = degree, raw = TRUE), data = data)
  } else if (ncol(X) == 4) {
    colnames(data) <- c("x1", "x2", "x3", "x4", "y")
    model <- lm(y ~ polym(x1, x2, x3, x4, degree = degree, raw = TRUE), data = data)
  } else if (ncol(X) == 5) {
    colnames(data) <- c("x1", "x2", "x3", "x4", "x5", "y")
    model <- lm(y ~ polym(x1, x2, x3, x4, x5, degree = degree, raw = TRUE), data = data)
  } else {
    stop("The function only takes 2 - 5 factors.")
  }
  # return
  model
}

# -------------------------------------------------------------------------

coef_name_to_array_index <- function(coefficients_name, k) {
  array_index_string <- stringr::str_extract(coefficients_name, "(\\d\\.)+[\\d]")
  array_index_number <- matrix(NA, length(array_index_string), k)
  array_index_number[1, ] <- 1
  for (i in 2:length(array_index_string)) {
    array_index_number[i, ] <- array_index_string[i] %>%
      stringr::str_split("\\.") %>%
      unlist() %>%
      as.numeric() + 1
  }
  # return
  array_index_number
}

# -------------------------------------------------------------------------

Gloptipolym <- function(coefficients, k, degree, lb, ub, maximization) {
  Ps <- list() # argument for GloptiPolyR, a list of lists
  # Objective function 
  P <- list()
  c <- array(0, dim = rep(degree + 1, k))
  # get position indices for the coefficients of the objective function
  id <- coef_name_to_array_index(names(coefficients), k = k)
  # put coefficient values into the multi-dimensional array
  for (i in 1:nrow(id)) {
    eval(parse(text = paste(
      "c[", toString(id[i, ]),
      "] <- coefficients[", i, "]"
    )))
    # example 1: eval(parse(text = "1+1")) -> 2
    # example 2: toString(id[1,]) -> "1, 1, 1"
  }
  if (maximization) { # assume GloptiPolyR only supports "min"
    P$c <- -c
  } else {
    P$c <- c
  }
  P$t <- "min" # specify attribute
  Ps[[1]] <- P # add to list
  # Constraint functions 
  for (i in 1:k) { # loop through variables
    # Lower bound constraint
    P <- list()
    c <- array(0, dim = rep(degree + 1, k))
    # specify coefficient 1 of the variable
    index_for_c <- rep(1, k)
    index_for_c[i] <- 2
    eval(parse(text = paste("c[", toString(index_for_c), "] <- 1")))
    # specify the constraint constant
    eval(parse(text = paste("c[", toString(rep(1, k)), "] <- -lb[", i, "]")))
    P$c <- c
    P$t <- ">=" # specify attribute
    Ps[[2 * i]] <- P # add to list
    # Upper bound constraint
    P <- list()
    c <- array(0, dim = rep(degree + 1, k))
    # specify coefficient 1 of the variable
    index_for_c <- rep(1, k)
    index_for_c[i] <- 2
    eval(parse(text = paste("c[", toString(index_for_c), "] <- 1")))
    # specify the constraint constant
    eval(parse(text = paste("c[", toString(rep(1, k)), "] <- -ub[", i, "]")))
    P$c <- c
    P$t <- "<="
    Ps[[2 * i + 1]] <- P
  }
  # Call GloptiPolyR 
  res <- GloptiPolyR(Ps)
  solution <- res$solution
  if (maximization) { # assume GloptiPolyR only supports "min"
    objective <- -res$objective
  } else {
    objective <- res$objective
  }
  # return
  list(solution = solution, objective = objective)
}
