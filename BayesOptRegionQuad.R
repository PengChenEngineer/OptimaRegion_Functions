BayesOptRegionQuad <- function(design, y, alpha = 0.05, n_post,
                               constr_lb, constr_ub, parallel = TRUE) {
  # draw posterior betas:
  c(betas_post, beta_hat) %<-% draw_post_quad_2(design = design, y = y, n_post = n_post)
  # trim posterior betas
  betas_post <- t(betas_post)
  d <- vector(length = nrow(betas_post))
  d <- DepthProc::depthTukey(betas_post, betas_post, ndir = 3000)
  order_d <- order(d)
  ind_alpha <- alpha * nrow(betas_post) + 1
  indices <- order_d[ind_alpha:nrow(betas_post)]
  betas_post <- betas_post[indices, ]
  betas_post <- t(betas_post)
  # map remaining posterior betas to optima
  x0s <- lhs::geneticLHS(n = 5, k = 2, criterium = "Maximin") %>%
    sweep(MARGIN = 2, STATS = constr_ub - constr_lb, FUN = "*") %>%
    sweep(MARGIN = 2, STATS = constr_lb, FUN = "+")
  if (parallel) {
    optima_post <- beta2opt_parApply(
      betas_post = betas_post, constr_lb = constr_lb, constr_ub = constr_ub, x0s = x0s
    )
  } else {
    optima_post <- beta2opt_apply(
      betas_post = betas_post, constr_lb = constr_lb, constr_ub = constr_ub, x0s = x0s
    )
  }
  optimum_post_mean <- beta2opt(
    beta = beta_hat, constr_lb = constr_lb, constr_ub = constr_ub, x0s = x0s
  )
  # return
  optima_post <- data.frame(optima_post)
  optimum_post_mean <- data.frame(optimum_post_mean)
  beta_hat <- t(beta_hat)
  beta_hat <- data.frame(beta_hat)
  beta_hat <- unlist(beta_hat[1, ])
  structure(
    list(
      optima = optima_post, opt_hat = optimum_post_mean,
      beta_hat = beta_hat,
      constr_lb = constr_lb, constr_ub = constr_ub
    ),
    class = "bayescrquad"
  )
}

# -------------------------------------------------------------------------

draw_post_quad_2 <- function(design, y, n_post) {
  X <- as.data.frame(design) %>%
    pmap_df(~ tibble(
      intercept = 1, x1 = ..1, x2 = ..2,
      x1x1 = ..1^2, x2x2 = ..2^2, x1x2 = ..1 * ..2
    )) %>%
    as.matrix(nrow = nrow(design)) # assume full quadratic basis
  beta_hat <- solve(crossprod(X), crossprod(X, y))
  s2 <- crossprod(y - X %*% beta_hat) / (nrow(design) - ncol(X))
  betas_post <- mvtnorm::rmvt(
    n = n_post,
    sigma = as.numeric(s2) * solve(crossprod(X)), df = nrow(design) - ncol(X),
    delta = beta_hat
  ) %>% t() # use argument delta instead of mu to avoid fallacy
  list(betas_post = betas_post, beta_hat = beta_hat)
}

# -------------------------------------------------------------------------

beta2opt_apply <- function(betas_post, constr_lb, constr_ub, x0s) {
  apply(betas_post, 2, beta2opt, constr_lb = constr_lb, constr_ub = constr_ub, x0s = x0s) %>%
    map_df(~ tibble(x1 = .x$x1, x2 = .x$x2)) %>%
    as.matrix()
}

# -------------------------------------------------------------------------

beta2opt_parApply <- function(betas_post, constr_lb, constr_ub, x0s) {
  cl <- makeCluster(
    detectCores(logical = FALSE) - 1,
    type = "FORK"
  )
  optima_post <- parApply(cl, betas_post, 2, beta2opt,
    constr_lb = constr_lb, constr_ub = constr_ub, x0s = x0s
  ) %>%
    map_df(~ tibble(x1 = .x$x1, x2 = .x$x2)) %>%
    as.matrix()
  stopCluster(cl)
  optima_post
}

# -------------------------------------------------------------------------

beta2opt <- function(beta, constr_lb, constr_ub, x0s) {
  nloptr_quad_2 <- purrr::partial(
    nloptr::nloptr,
    eval_f = function(x) -quad_2(x, beta),
    eval_grad_f = function(x) -quad_2_grad(x, beta),
    lb = constr_lb, ub = constr_ub,
    opts = list("algorithm" = "NLOPT_LD_MMA", print_level = 0, xtol_rel = 1e-03)
  )
  lapply(1:nrow(x0s), function(i) x0s[i, ]) %>%
    purrr::map(function(i) nloptr_quad_2(x0 = i)) %>%
    map_df(~ tibble(
      x1 = .x$solution[1], x2 = .x$solution[2], objective = -.x$objective
    )) %>%
    top_n(1, objective) %>%
    distinct() %>%
    "["(1, ) %>%
    dplyr::select(x1, x2)
}

# -------------------------------------------------------------------------

quad_2 <- function(x, beta) {
  beta[1] + beta[2] * x[1] + beta[3] * x[2] +
    beta[4] * x[1]^2 + beta[5] * x[2]^2 + beta[6] * x[1] * x[2]
}

# -------------------------------------------------------------------------

quad_2_grad <- function(x, beta) {
  rbind(
    beta[2] + 2 * beta[4] * x[1] + beta[6] * x[2],
    beta[3] + 2 * beta[5] * x[2] + beta[6] * x[1]
  )
}

# -------------------------------------------------------------------------

plot.bayescrquad <- function(res, xlab = "x1", ylab = "x2") {
  # contour grid
  grid_size <- 200
  x1_grid <- seq(res$constr_lb[1], res$constr_ub[1], length.out = grid_size)
  x2_grid <- seq(res$constr_lb[2], res$constr_ub[2], length.out = grid_size)
  contour_grid <- expand.grid(x1_grid, x1_grid)
  y_grid <- apply(contour_grid, 1, quad_2, res$beta_hat)
  contour_data <- data.frame(
    x1 = contour_grid[, 1], x2 = contour_grid[, 2], y = y_grid
  )
  # plot
  ggplot(data = contour_data, aes(x = x1, y = x2)) +
    geom_polygon(
      data = res$optima %>% slice(chull(x1, x2)),
      mapping = aes(x = x1, y = x2),
      alpha = 0.8, col = "black", size = 0.3, fill = "gray"
    ) +
    geom_contour(aes(z = y), col = "black", size = 0.5) +
    geom_text_contour(aes(z = y), stroke = 0.1) +
    geom_point(data = res$opt_hat, aes(x = x1, y = x2), col = "red") +
    coord_fixed() +
    xlab(xlab) + ylab(ylab) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}