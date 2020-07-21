OptRegionPC1Quad <- function(design, Y, constr_lb, constr_ub, alpha = 0.05, B = 200) {
  X <- design
  ppca <- my_ppca(X = Y, M = 1)
  cr <- OptRegionQuad(
    design = X, y = ppca$z_given_x[1, ], constr_lb = constr_lb, constr_ub = constr_ub,
    alpha = alpha, B = B
  )
  list(cr = cr, ppca = ppca)
}


# -------------------------------------------------------------------------

my_ppca <- function(X, M) {
  D = nrow(X)
  N = ncol(X)
  # eigen operation
  x_bar <- apply(X, 1, mean) # (D, 1)
  X_center <- sweep(X, 1, x_bar) # (D, N)
  S <- tcrossprod(X_center) / N # (D, D)
  S_eigen <- eigen(S)
  props <-  S_eigen$values / sum(S_eigen$values)
  if(sum(S_eigen$vectors[, 1]) < 0) 
    S_eigen$vectors[, 1] <- - S_eigen$vectors[, 1] # adjust pc1 direction
  # mle for mu
  mu_hat <- x_bar 
  # mle for sigma2
  if(M == D) {
    sigma2_hat <- 0
  } else {
    sigma2_hat <- sum(S_eigen$values[-(1:M)]) / (D - M)
  }
  # mle for W
  U <- S_eigen$vectors[, 1:M] %>% matrix(nrow = D) # (D, M)
  if(M == 1) {
    other_term <- sqrt(S_eigen$values[1] - sigma2_hat) # (1, 1)
    W_hat <- U * other_term # (D, 1)
  } else {
    other_term <- (S_eigen$values[1:M] - sigma2_hat) %>% sqrt() %>% diag() # (M, M)
    W_hat <- U %*% other_term # (D, M)
  }
  # z|x
  M_hat <- crossprod(W_hat) + sigma2_hat * diag(M) # (M, M) 
  z_given_x <- solve(M_hat, crossprod(W_hat, X_center)) # (M, N)
  # return
  list(
    mu_hat = mu_hat, sigma2_hat = sigma2_hat, W_hat = W_hat,
    z_given_x = z_given_x, props = round(props, 4)
  )
}

