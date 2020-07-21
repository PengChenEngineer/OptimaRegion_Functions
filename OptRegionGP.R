OptRegionGP <- function(design, y, constr_lb, constr_ub, alpha = 0.05, B = 200, 
                         xi, base_grid_size = 225, parallel = TRUE){
  # rs_hat and opt_hat
  message("Finding point estimates of the response surface and its optimum ...")
  mod_km_hat <- DiceKriging::km(
    design = design, response = y, covtype = "exp",
    iso = TRUE, nugget.estim = TRUE, control = list(trace = FALSE)
  )
  theta_hat <- mod_km_to_theta(mod_km_hat)
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
        sample_path_smoother(
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
        sample_path_smoother(
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

plot.crgpok <- function(res, xlab = "x1", ylab = "x2") {
  optima <- res$optima
  opt_bag <- res$opt_bag
  x2krig_mean <- res$rs_hat
  constr_lb <- res$constr_lb
  constr_ub <- res$constr_ub
  # contour grid
  grid_size <- 100
  x1_grid <- seq(constr_lb[1], constr_ub[1], length.out = grid_size)
  x2_grid <- seq(constr_lb[2], constr_ub[2], length.out = grid_size)
  contour_grid <- expand.grid(x1_grid, x1_grid)
  y_grid <- apply(contour_grid, 1, x2krig_mean)
  contour_data <- data.frame(
    x1 = contour_grid[, 1], x2 = contour_grid[, 2], y = y_grid
  )
  ggplot(data = contour_data, aes(x = x1, y = x2)) + 
    geom_polygon(
      data = optima %>% slice(chull(x1, x2)),
      mapping = aes(x = x1, y = x2),
      alpha = 0.8, col = "black", size = 0.3, fill = "gray"
    ) + 
    geom_contour(aes(z = y), col = "black", size = 0.5) +
    geom_text_contour(aes(z = y), stroke = 0.1) +
    geom_point(
      data = tibble(x = opt_bag$x1, y = opt_bag$x2),
      mapping = aes(x, y),
      col = "red"
    )+
    coord_fixed() +
    xlab(xlab) + ylab(ylab) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    )
}

# -------------------------------------------------------------------------

mod_km_to_theta <- function(mod_km) {
  list(
    process_mean = mod_km@trend.coef,
    length_scale = mod_km@covariance@range.val, 
    fun_var = mod_km@covariance@sd2,
    noise_var = mod_km@covariance@nugget
  )
}

# -------------------------------------------------------------------------

theta_to_obs <- function(theta, points, y) {
  kernel <- function(d) exp(-d / theta$length_scale)
  R11 <- apply(X = fields::rdist(points), MARGIN = c(1, 2), FUN = kernel)
  mu <- rep(theta$process_mean, nrow(points))
  Sigma <- theta$fun_var * R11 + theta$noise_var * diag(nrow(R11))
  MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma)
}

# -------------------------------------------------------------------------

sample_path_smoother <- function(theta_hat, design, y, all_points,
                                 decomp_1st_term, x0s, 
                                 constr_lb, constr_ub) {
  ## noisy observation at old points
  obs_star <- theta_to_obs(
    theta = theta_hat, points = design, y = y
  )
  ## latent process at all points
  process_star <- theta2paths(
    theta = theta_hat, n_path_per_theta = 1, path_grid = all_points,
    design = design, y = obs_star
  ) %>% `[[`(1)
  ## smoother from old points at all points
  mod_km_star_1 <- DiceKriging::km(
    design = design, response = obs_star, covtype = "exp",
    iso = TRUE, nugget.estim = TRUE, control = list(trace = FALSE)
  )
  theta_hat_star_1 <- mod_km_to_theta(mod_km_star_1)
  rs_hat_star_1 <- param2krig_fn(
    design = design, y = obs_star, theta = theta_hat_star_1
  )
  smoother_vec_1 <- purrr::map_dbl(
    .x = 1:nrow(all_points), .f = ~ rs_hat_star_1(all_points[.x, ])
  )
  ## simulated second term of decomposition
  decomp_2nd_term <- process_star - smoother_vec_1
  ## sample-path smoother and its optimum
  sample_path_vec <- decomp_1st_term + decomp_2nd_term
  mod_km_sample_path <- DiceKriging::km(
    design = all_points, response = sample_path_vec, covtype = "exp",
    iso = TRUE, nugget.estim = TRUE, control = list(trace = FALSE)
  )
  theta_hat_sample_path <- mod_km_to_theta(mod_km_sample_path)
  sample_path_smoother <- param2krig_fn(
    design = all_points, y = sample_path_vec, theta = theta_hat_sample_path
  )
  sample_path_optimum <- krig_fn2optimum(
    krig_fn = sample_path_smoother,
    constr_lb = constr_lb, constr_ub = constr_ub, x0s = x0s
  )
  ## return
  list(
    theta_hat_sample_path = theta_hat_sample_path %>% unlist(),
    sample_path_optimum = sample_path_optimum %>% unlist()
  )
}
