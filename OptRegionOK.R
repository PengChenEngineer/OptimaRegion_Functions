OptRegionOK <- function(design, y, constr_lb, constr_ub, alpha = 0.05, B = 200, 
                         xi, base_grid_size = 225, parallel = TRUE) {
  # rs_hat and opt_hat
  message("Finding point estimates of the response surface and its optimum ...")
  theta_hat <- fit_variogram(design = design, y = y)
  rs_hat <- param2krig_fn(design = design, y = y, theta = theta_hat)
  x0s <- lhs::geneticLHS(n = 20, k = ncol(design), criterium = "Maximin") %>%
    sweep(MARGIN = 2, STATS = constr_ub - constr_lb, FUN = "*") %>%
    sweep(MARGIN = 2, STATS = constr_lb, FUN = "+")
  opt_hat <- krig_fn2optimum(
    krig_fn = rs_hat, constr_lb = constr_lb, constr_ub = constr_ub, x0s = x0s
  )
  # bootstrap
  message("Bootstrapping ...")
  ## old and new points
  all_points <- generate_path_grid(
    xi = xi, base_grid_size = base_grid_size, design = design, opt_hat = opt_hat,
    constr_lb = constr_lb, constr_ub = constr_ub
  )
  ## first term of decomposition
  decomp_1st_term <- purrr::map_dbl(
    .x = 1:nrow(all_points), .f = ~ rs_hat(all_points[.x, ])
  )
  # map to B approximated sample paths
  if(parallel) {
    cl <- makeCluster(
      detectCores(logical = FALSE) - 1, outfile = "log.txt", type = "FORK"
    )
    sample_paths_approx <- clusterApply(
      cl = cl,
      x = 1:B,
      fun = function(i) {
        if(i %% 10 == 0) paste0("The ", i, " th bootstrap ...") %>% print()
        sample_path_smoother_OK(
          theta_hat = theta_hat, design = design, y = y, all_points = all_points,
          decomp_1st_term = decomp_1st_term, x0s = x0s,
          constr_lb = constr_lb, constr_ub = constr_ub
        )
      }
    )
    stopCluster(cl)
  } else {
    sample_paths_approx <- lapply(
      X = 1:B,
      FUN = function(i) {
        if(i %% 10 == 0) paste0("The ", i, " th bootstrap ...") %>% print()
        sample_path_smoother_OK(
          theta_hat = theta_hat, design = design, y = y, all_points = all_points,
          decomp_1st_term = decomp_1st_term, x0s = x0s,
          constr_lb = constr_lb, constr_ub = constr_ub
        )
      }
    )
  }
  # matrices for thetas and optima
  thetas <- matrix(NA, nrow = B, ncol = 4)
  optima <- matrix(NA, nrow = B, ncol = 2)
  for(b in 1:B) {
    thetas[b, ] <- sample_paths_approx[[b]]$theta_hat_sample_path
    optima[b, ] <- sample_paths_approx[[b]]$sample_path_optimum
  }
  # obtain final optima
  message("Computing the confidence region ...")
  ## trim optima according to thetas
  d <- DepthProc::depthTukey(u = thetas, X = thetas, ndir = 3000)
  order_d <- order(d)
  indices <- order_d[(alpha * B + 1):B]
  optima <- optima[indices, ]
  ## remove optima outside convex hull
  optima_inCR <- purrr::map_lgl(
    .x = 1:nrow(optima),
    .f = ~ inCR(new_point = optima[.x, ], old_points = design)
  )
  optima <- optima[optima_inCR, ]
  # return 
  optima <- data.frame(x1 = optima[, 1], x2 = optima[, 2])
  opt_bag <- apply(optima, 2, mean) 
  opt_bag <- data.frame(x1 = opt_bag[1], x2 = opt_bag[2])
  structure(
    list(
      optima = optima, opt_bag = opt_bag,
      rs_hat = rs_hat,
      constr_lb = constr_lb, constr_ub = constr_ub
    ),
    class = "crgpok"
  )
}

# -------------------------------------------------------------------------

fit_variogram <- function(design, y) {
  design_matrix <- design
  responses <- y
  data <- data.frame(design_matrix, responses)
  colnames(data) <- c("x1", "x2", "y")
  sp::coordinates(data) <- ~ x1 + x2
  empirical_variogram <- gstat::variogram(y ~ 1, data = data)
  variogram_mod <- gstat::vgm(
    model = "Exp",
    range = (max(design_matrix[, 1]) - min(design_matrix[, 1]) +
      max(design_matrix[, 2]) - min(design_matrix[, 2])) / 20,
    psill = (max(responses) - min(responses)) / 10,
    nugget = (max(responses) - min(responses)) / 100
  )
  suppressWarnings(
    exp_variogram <- gstat::fit.variogram(empirical_variogram, variogram_mod)
  )
  # extract cov params
  length_scale = exp_variogram$range[2]
  fun_var = exp_variogram$psill[2]
  noise_var = exp_variogram$psill[1]
  # ones
  ones <- matrix(rep(1, nrow(design)), nrow = nrow(design))
  # precision matrix
  kernel <- function(d) exp(-d / length_scale)
  R11 <- apply(X = fields::rdist(design), MARGIN = c(1, 2), FUN = kernel)
  Sigma <- fun_var * R11 + noise_var * diag(nrow(R11))
  precise_mat <- solve(Sigma)
  # process mean
  process_mean_n <- crossprod(ones, precise_mat) %*% y
  process_mean_d <- crossprod(ones, precise_mat) %*% ones
  process_mean <- as.numeric(process_mean_n / process_mean_d)
  # return
  list(
    process_mean = process_mean, length_scale = length_scale,
    fun_var = fun_var, noise_var = noise_var,
    ones = ones, precise_mat = precise_mat, kernel = kernel
  )
}

# -------------------------------------------------------------------------

sample_path_smoother_OK <- function(theta_hat, design, y, all_points,
                                    decomp_1st_term, x0s, 
                                    constr_lb, constr_ub) {
  ## process at all points
  ### decor and center y
  center_mat <- diag(nrow(design)) - tcrossprod(theta_hat$ones) / nrow(design)
  z_X <- center_mat %*% expm::sqrtm(theta_hat$precise_mat) %*% y 
  z_X <- as.numeric(z_X)
  ### resample
  z_XS_star <- sample(x = z_X, size = nrow(all_points), replace = TRUE)
  ### recor and add mean
  K_XS_XS <- apply(
    X = fields::rdist(all_points), MARGIN = c(1, 2), FUN = theta_hat$kernel
  )
  f_XS_star <- matrix(rep(theta_hat$process_mean, nrow(all_points)), ncol = 1) + 
    sqrt(theta_hat$fun_var) * expm::sqrtm(K_XS_XS) %*% z_XS_star
  f_XS_star <- as.numeric(f_XS_star)
  ## observations at old points
  y_X_star <- f_XS_star[1:nrow(design)] +
    sqrt(theta_hat$noise_var) * z_XS_star[1:nrow(design)]
  ## smoother from old points at all points
  theta_hat_star_1 <- fit_variogram(design = design, y = y_X_star)
  rs_hat_star_1 <- param2krig_fn(
    design = design, y = y_X_star, theta = theta_hat_star_1
  )
  smoother_vec_1 <- purrr::map_dbl(
    .x = 1:nrow(all_points), .f = ~ rs_hat_star_1(all_points[.x, ])
  )
  ## simulated second term of decomposition
  decomp_2nd_term <- f_XS_star - smoother_vec_1
  ## sample-path smoother and its optimum
  sample_path_vec <- decomp_1st_term + decomp_2nd_term
  theta_hat_sample_path <- fit_variogram(design = all_points, y = sample_path_vec)
  sample_path_smoother <- param2krig_fn(
    design = all_points, y = sample_path_vec, theta = theta_hat_sample_path
  )
  sample_path_optimum <- krig_fn2optimum(
    krig_fn = sample_path_smoother,
    constr_lb = constr_lb, constr_ub = constr_ub, x0s = x0s
  )
  # return
  list(
    theta_hat_sample_path = theta_hat_sample_path[1:4] %>% unlist(),
    sample_path_optimum = sample_path_optimum %>% unlist()
  )
}
