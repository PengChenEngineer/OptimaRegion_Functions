BayesOptRegionGP <- function(design, y, constr_lb, constr_ub, alpha = 0.05,
                              chain_length, thin_interval, n_path_per_theta,
                              xi, base_grid_size = 225,
                              process_mean_lb, process_mean_ub,
                              length_scale_lb, length_scale_ub,
                              fun_var_lb, fun_var_ub, noise_var_lb, noise_var_ub,
                              parallel = TRUE) {
  # draw thetas
  message("simulating posterior parameters ... ")
  thetas <- draw_post_gp(
    design = design, y = y,
    process_mean_lb = process_mean_lb, process_mean_ub = process_mean_ub,
    length_scale_lb = length_scale_lb, length_scale_ub = length_scale_ub,
    fun_var_lb = fun_var_lb, fun_var_ub = fun_var_ub,
    noise_var_lb = noise_var_lb, noise_var_ub = noise_var_ub,
    chain_length = chain_length
  )
  # generate sample paths grid
  theta_mean <- purrr::map(thetas, mean)
  x2krig_mean <- param2krig_fn(design = design, y = y, theta = theta_mean)
  opt_hat <- krig_fn2optimum(
    krig_fn = x2krig_mean, constr_lb = constr_lb, constr_ub = constr_ub, x0s = design
  )
  path_grid <- generate_path_grid(
    xi = xi, base_grid_size = base_grid_size, design = design, opt_hat = opt_hat,
    constr_lb = constr_lb, constr_ub = constr_ub
  )
  # thin thetas
  thetas_thin <- slice(thetas, seq(0, nrow(thetas), by = thin_interval))
  # map thetas to smoothers
  message("mapping posterior parameters to smoothers ... ") 
  smoothers <- thetas2smoothers(
    thetas = thetas_thin, n_path_per_theta = n_path_per_theta, 
    path_grid = path_grid, design = design, y = y
  ) %>%
    unlist()
  # trim smoothers
  message("trimming the smoothers ... ")
  alphas <- smoothers2alphas(smoothers = smoothers) 
  alphas <- do.call(rbind, alphas)
  d <- vector(length = nrow(alphas))
  d <- DepthProc::depthTukey(alphas, alphas, ndir = 3000)
  order_d <- order(d)
  ind_alpha <- alpha * nrow(alphas) + 1
  indices <- order_d[ind_alpha:nrow(alphas)]
  smoothers <- smoothers[indices]
  # map smoothers to optima
  message("mapping smoothers to optima ... ")
  if (parallel) {
    optima <- smoothers2optima_par(
      smoothers = smoothers, design = design,
      constr_lb = constr_lb, constr_ub = constr_ub
    )
  } else {
    optima <- smoothers2optima(
      smoothers = smoothers, design = design, 
      constr_lb = constr_lb, constr_ub = constr_ub
    )
  }
  # return
  opt_hat <- data.frame(opt_hat)
  optima <- data.frame(optima)
  list(
    thetas = structure(thetas, class = "mcmcdraw"),
    thetas_thin = structure(thetas_thin, class = "mcmcdraw"),
    credible_region = structure(
      list(
        optima = optima,
        opt_hat = opt_hat,
        rs_hat = x2krig_mean,
        constr_lb = constr_lb,
        constr_ub = constr_ub
      ),
      class = "bayescrgp"
    )
  )
}

# -------------------------------------------------------------------------

draw_post_gp <- function(design, y, process_mean_lb, process_mean_ub,
                         length_scale_lb, length_scale_ub,
                         fun_var_lb, fun_var_ub, noise_var_lb, noise_var_ub,
                         chain_length) {
  # define the bayesian model
  gp_model <- "model{
    # likelihood
    y ~ dmnorm(mu, inverse(Sigma))
      for(i in 1:length(y)) {
        mu[i] <- process_mean
      }
      for(i in 1:(length(y) - 1)) {
        for(j in (i + 1):length(y)) {
          Sigma[i, j] <- fun_var * exp(- d[i, j] / length_scale)
          Sigma[j, i] <- Sigma[i, j]
        }
      }
      for(i in 1:length(y)) {
        Sigma[i, i] <- fun_var + noise_var
      }
    # prior
    process_mean ~ dunif(process_mean_lb, process_mean_ub)
    length_scale ~ dunif(length_scale_lb, length_scale_ub)
    fun_var ~ dunif(fun_var_lb, fun_var_ub)
    noise_var ~ dunif(noise_var_lb, noise_var_ub)
  }"
  # compile the bayesian model
  gp_jags <- jags.model(
    textConnection(gp_model),
    data = list(
      y = y, d = fields::rdist(design),
      process_mean_lb = process_mean_lb, process_mean_ub = process_mean_ub,
      length_scale_lb = length_scale_lb, length_scale_ub = length_scale_ub,
      fun_var_lb = fun_var_lb, fun_var_ub = fun_var_ub,
      noise_var_lb = noise_var_lb, noise_var_ub = noise_var_ub
    ),
    inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100),
    quiet = TRUE
  )
  # simulate posterior parameters
  gp_sim <- coda.samples(
    model = gp_jags,
    variable.names = c("process_mean", "length_scale", "fun_var", "noise_var"),
    n.iter = chain_length
  )
  # extract useful 
  gp_chain <- data.frame(gp_sim[[1]]) %>%
    dplyr::select(process_mean, length_scale, fun_var, noise_var)
  # return
  return(gp_chain)
}

# -------------------------------------------------------------------------

param2krig_fn <- function(design, y, theta) {
  kernel <- function(d) exp(-d / theta$length_scale)
  R11 <- apply(
    X = fields::rdist(design), MARGIN = c(1, 2), FUN = kernel
  )
  alpha <- theta$fun_var * solve(
    theta$fun_var * R11 + theta$noise_var * diag(nrow(R11)),
    y - rep(theta$process_mean, length(y))
  ) 
  rm(R11)
  
  function(x) {
    x <- matrix(x, nrow = 1)
    R12 <- apply(
      X = fields::rdist(x1 = design, x2 = x), MARGIN = c(1, 2), FUN = kernel
    )
    theta$process_mean + crossprod(R12, alpha)
  }
}

# -------------------------------------------------------------------------

krig_fn2optimum <- function(krig_fn, constr_lb, constr_ub, x0s) {
  x0s <- as.matrix(x0s)
  nloptr_gp <- purrr::partial(
    nloptr::nloptr,
    eval_f = function(x) -krig_fn(x),
    lb = constr_lb, ub = constr_ub,
    opts = list("algorithm" = "NLOPT_LN_COBYLA", print_level = 0, "xtol_rel" = 1.0e-3)
  )
  
  purrr::map(
    .x = 1:nrow(x0s),
    .f = ~ nloptr_gp(x0s[.x, ])
  ) %>%
    purrr::map_df(
      .f = ~ tibble(x1 = .x$solution[1], x2 = .x$solution[2], objective = -.x$objective)
    ) %>%
    top_n(1, objective) %>%
    distinct() %>%
    dplyr::select(x1, x2)
}

# -------------------------------------------------------------------------

generate_path_grid <- function(xi, base_grid_size, design, opt_hat,
                               constr_lb, constr_ub) {
  base_grid <- expand.grid(
    x1 = seq(constr_lb[1], constr_ub[1], length.out = sqrt(base_grid_size)),
    x2 = seq(constr_lb[2], constr_ub[2], length.out = sqrt(base_grid_size))
  ) %>% as.matrix()
  gird_inCR <- purrr::map_lgl(
    .x = 1:nrow(base_grid),
    .f = ~ inCR(new_point = base_grid[.x, ], old_points = design)
  )
  candidate_grid <- base_grid[gird_inCR, ] 
  
  opt_hat <- matrix(opt_hat, nrow = 1)
  path_grid <- design # include old points 
  n_new <- nrow(candidate_grid) / 5 %>% as.integer() # implicit M
  for (i in 1:n_new) {
    scores <- purrr::map_dbl(
      .x = 1:nrow(candidate_grid),
      .f = ~ min(
        fields::rdist(
          x1 = candidate_grid[.x, ] %>% matrix(nrow = 1), x2 = path_grid
        )
      ) -
        xi * fields::rdist(
          x1 = candidate_grid[.x, ] %>% matrix(nrow = 1), x2 = opt_hat
        )
    )
    s_star <- candidate_grid[which.max(scores), ]
    candidate_grid <- candidate_grid[-which.max(scores), ]
    path_grid <- rbind(path_grid, s_star) 
  }
  path_grid
}

# -------------------------------------------------------------------------

inCR <- function(new_point, old_points) {
  old_points <- as.matrix(old_points)
  dim(new_point) <- c(1, 2)
  # indices of the points lying the convex hull
  coords <- grDevices::chull(old_points)
  # area of the convex hull
  CR_area <- geometry::polyarea(old_points[coords, 1], old_points[coords, 2])
  # add the new point to the set and re-calculate the area
  r <- rbind(old_points, new_point)
  newCH_coords <- chull(r)
  new_area <- geometry::polyarea(r[newCH_coords, 1], r[newCH_coords, 2])
  # tell if the new point is in the original convex hull by comparing the areas
  if (new_area > CR_area) FALSE else TRUE
}

# -------------------------------------------------------------------------

thetas2smoothers <- function(thetas, n_path_per_theta, path_grid, design, y) {
  purrr::map(
    .x = 1:nrow(thetas),
    .f = ~ theta2smoothers(
      theta = thetas[.x, ], n_path_per_theta = n_path_per_theta, path_grid = path_grid,
      design = design, y = y
    )
  )
}

# -------------------------------------------------------------------------

theta2smoothers <- function(theta, n_path_per_theta, path_grid, design, y) {
  paths <- theta2paths(
    theta = theta, n_path_per_theta = n_path_per_theta, path_grid = path_grid,
    design = design, y = y
  )
  x2krigs <- purrr::map(
    .x = paths,
    .f = ~ param2krig_fn(
      design = path_grid, y = .x, theta = theta
    )
  )
  x2krigs
}

# -------------------------------------------------------------------------

theta2paths <- function(theta, n_path_per_theta, path_grid, design, y) {
  kernel <- function(d) exp(-d / theta$length_scale)
  R11 <- apply(X = fields::rdist(design), MARGIN = c(1, 2), FUN = kernel)
  R12 <- apply(
    X = fields::rdist(x1 = design, x2 = path_grid),
    MARGIN = c(1, 2), FUN = kernel
  )
  R22 <- apply(X = fields::rdist(path_grid), MARGIN = c(1, 2), FUN = kernel)
  mu <- rep(theta$process_mean, nrow(path_grid)) +
    crossprod(
      theta$fun_var * R12,
      solve(
        theta$fun_var * R11 + theta$noise_var * diag(nrow(R11)),
        y - rep(theta$process_mean, length(y))
      )
    )
  Sigma <- theta$fun_var * R22 -
    crossprod(
      theta$fun_var * R12,
      solve(
        theta$fun_var * R11 + theta$noise_var * diag(nrow(R11)),
        theta$fun_var * R12
      )
    )
  MASS::mvrnorm(n = n_path_per_theta, mu = mu, Sigma = Sigma) %>%
    matrix(nrow = n_path_per_theta) %>%
    t() %>%
    as.data.frame() %>%
    set_names(map_chr(1:n_path_per_theta, ~ paste("path_", .x, sep = "")))
}

# -------------------------------------------------------------------------

smoothers2alphas <- function(smoothers) {
  purrr::map(
    .x = smoothers,
    .f = ~ smoother2alpha(.x)
  )
}

# -------------------------------------------------------------------------

smoother2alpha <- function(smoother){
  smoother_envir <- environment(smoother)
  smoother_envir$alpha
}

# -------------------------------------------------------------------------

smoothers2optima <- function(smoothers, design, constr_lb, constr_ub) {
  x0s <- lhs::geneticLHS(n = 20, k = ncol(design), criterium = "Maximin") %>%
    sweep(MARGIN = 2, STATS = constr_ub - constr_lb, FUN = "*") %>%
    sweep(MARGIN = 2, STATS = constr_lb, FUN = "+") %>%
    `colnames<-`(c("x1", "x2"))
  optima <- lapply(
    X = 1:length(smoothers),
    FUN = function(i) {
      paste0("optimizing the ", i, " th smoother ...") %>%
        print()
      krig_fn2optimum(
        krig_fn = smoothers[[i]], constr_lb = constr_lb, constr_ub = constr_ub, x0s = x0s
      )
    }
  )
  optima <- purrr::map_df(
    .x = optima,
    .f = ~ tibble(x1 = .x$x1, x2 = .x$x2)
  ) %>%
    as.matrix()
  optima_inCR <- purrr::map_lgl(
    .x = 1:nrow(optima),
    .f = ~ inCR(new_point = optima[.x, ], old_points = design)
  )
  optima[optima_inCR, ]
}

# -------------------------------------------------------------------------

smoothers2optima_par <- function(smoothers, design, constr_lb, constr_ub) {
  x0s <- lhs::geneticLHS(n = 20, k = ncol(design), criterium = "Maximin") %>%
    sweep(MARGIN = 2, STATS = constr_ub - constr_lb, FUN = "*") %>%
    sweep(MARGIN = 2, STATS = constr_lb, FUN = "+") %>%
    `colnames<-`(c("x1", "x2"))
  # parallel starts
  cl <- makeCluster(
    detectCores(logical = FALSE) - 1, outfile = "log.txt", type = "FORK"
  )
  optima <- clusterApply(
    cl = cl,
    x = 1:length(smoothers),
    fun = function(i) {
      paste0("optimizing ", i, " th smoother ...") %>%
        print()
      krig_fn2optimum(
        krig_fn = smoothers[[i]], constr_lb = constr_lb, constr_ub = constr_ub, x0s = x0s
      )
    }
  )
  stopCluster(cl)
  # parallel ends
  optima <- purrr::map_df(
    .x = optima,
    .f = ~ tibble(x1 = .x$x1, x2 = .x$x2)
  ) %>%
    as.matrix()
  optima_inCR <- purrr::map_lgl(
    .x = 1:nrow(optima),
    .f = ~ inCR(new_point = optima[.x, ], old_points = design)
  )
  optima[optima_inCR, ]
}

# -------------------------------------------------------------------------

plot.bayescrgp <- function(credible_region, xlab = "x1", ylab = "x2") {
  class(credible_region) <- "list"
  c(optima, opt_hat, x2krig_mean, constr_lb, constr_ub) %<-% credible_region
  # contour grid
  grid_size <- 200
  x1_grid <- seq(constr_lb[1], constr_ub[1], length.out = grid_size)
  x2_grid <- seq(constr_lb[2], constr_ub[2], length.out = grid_size)
  contour_grid <- expand.grid(x1_grid, x1_grid)
  y_grid <- apply(contour_grid, 1, x2krig_mean)
  contour_data <- data.frame(
    x1 = contour_grid[, 1], x2 = contour_grid[, 2], y = y_grid
  )
  # extract design
  design <- fn_env(x2krig_mean)$design
  # plot
  ggplot(data = contour_data,aes(x = x1, y = x2)) +
    geom_polygon(
      data = optima %>% slice(chull(x1, x2)),
      mapping = aes(x = x1, y = x2),
      alpha = 0.8, col = "black", size = 0.3, fill = "gray"
    ) +
    geom_contour(aes(z = y), col = "black", size = 0.5) +
    geom_text_contour(aes(z = y), stroke = 0.1) +
    geom_point(
      data = tibble(x = opt_hat$x1, y = opt_hat$x2),
      mapping = aes(x, y),
      col = "red"
    ) +
    coord_fixed() +
    xlab(xlab) + ylab(ylab) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

# -------------------------------------------------------------------------

plot.mcmcdraw <- function(thetas) {
  attr(thetas, "class") <- "data.frame"
  traces <- purrr::imap(
    thetas,
    ~ ggplot(thetas, aes(x = 1:length(.x), y = .x)) +
      geom_line() +
      geom_smooth(method = "lm", col = "red") +
      xlab("Iterations") +
      ylab("") +
      ggtitle(paste("trace of posterior", .y)) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )
  )
  densities <- purrr::imap(
    thetas,
    ~ ggplot(thetas, aes(x = .x)) +
      geom_density() +
      xlab("Iterations") +
      ylab("") +
      ggtitle(paste("density of posterior", .y)) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )
  )
  grid.arrange(
    grobs = c(traces, densities),
    layout_matrix = matrix(1:8, nrow = 4)
  )
}

# -------------------------------------------------------------------------

summary.mcmcdraw <- function(thetas) {
  attr(thetas, "class") <- "data.frame"
  purrr::map(thetas, ~ quantile(x = .x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))) %>%
    unlist() %>%
    matrix(nrow = 4, byrow = TRUE) %>%
    `rownames<-`(names(thetas)) %>%
    `colnames<-`(c("2.5%", "25%", "50%", "75%", "97.5%"))
}