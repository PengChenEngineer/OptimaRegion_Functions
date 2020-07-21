#' Global optimization of up to cubic polynomial functions (up to 5 variables)
#'
#' Optimize a quadratic or cubic polynomial functionin 2 ~ 5 variables
#' with bound constraints
#' \insertCite{DelCastilloCR}{OptimaRegion}.
#'
#' GloptipolyR is an R implementation of the “Gloptipoly” algorithm
#' \insertCite{lasserre2001global}{OptimaRegion}
#'
#' @param P A list of list; Each sub-list has 2 components:
#' 1. a multi-dimensional array corresponding to a objective or constraint function
#' 2. an attribute of the objective or constraint function
#'
#' @return Returns the optimal solution and its corresponding objective value
#' @inheritSection OptRegionQuad Author(s)
#' @references{
#'  \insertAllCited{}
#' }
#' @examples
#' # Optimize the following quadratic function in 3 variables
#' # f(x) = -1.5 x_1 + 2.13 x_2 - 1.81 x_3 + 7.13 x_1 x_2 +
#' #        3.27 x_1 x_3  + 2.73 x_2 x_3 +
#' #        4.69 x_1^2 + 6.27 x_2^2 + 5.21 x_3^2.
#'
#' # The input for GloptiPolyR is a list of 7 sub-lists,
#' # each of which corresponds to the objective function or a constraint
#' # function, respectively. See del Castillo et al. (2019) for details.
#' P <- list()
#' p_f <- list()
#' p_g_1 <- list()
#' p_g_2 <- list()
#' p_g_3 <- list()
#' p_g_4 <- list()
#' p_g_5 <- list()
#' p_g_6 <- list()
#'
#' p_f$c <- array(0, dim = c(3, 3, 3))
#' p_f$c[2, 1, 1] <- -1.5
#' p_f$c[1, 2, 1] <- 2.13
#' p_f$c[1, 1, 2] <- -1.81
#' p_f$c[2, 2, 1] <- 7.13
#' p_f$c[2, 1, 2] <- 3.27
#' p_f$c[1, 2, 2] <- 2.73
#' p_f$c[3, 1, 1] <- 4.69
#' p_f$c[1, 3, 1] <- 6.27
#' p_f$c[1, 1, 3] <- 5.21
#'
#' p_g_1$c <- array(0, dim = c(3, 3, 3))
#' p_g_1$c[1, 1, 1] <- 2
#' p_g_1$c[2, 1, 1] <- 1
#'
#' p_g_2$c <- array(0, dim = c(3, 3, 3))
#' p_g_2$c[1, 1, 1] <- -2
#' p_g_2$c[2, 1, 1] <- 1
#'
#' p_g_3$c <- array(0, dim = c(3, 3, 3))
#' p_g_3$c[1, 1, 1] <- 2
#' p_g_3$c[1, 2, 1] <- 1
#'
#' p_g_4$c <- array(0, dim = c(3, 3, 3))
#' p_g_4$c[1, 1, 1] <- -2
#' p_g_4$c[1, 2, 1] <- 1
#'
#' p_g_5$c <- array(0, dim = c(3, 3, 3))
#' p_g_5$c[1, 1, 1] <- 2
#' p_g_5$c[1, 1, 2] <- 1
#'
#' p_g_6$c <- array(0, dim = c(3, 3, 3))
#' p_g_6$c[1, 1, 1] <- -2
#' p_g_6$c[1, 1, 2] <- 1
#'
#' # Set the attribute for the objective function as either ``min'' or ``max''.
#' p_f$t <- "min"
#'
#' # Set the attributes for the constraint functions as either ``>='' or ``<=''.
#' p_g_1$t <- ">="
#' p_g_2$t <- "<="
#' p_g_3$t <- ">="
#' p_g_4$t <- "<="
#' p_g_5$t <- ">="
#' p_g_6$t <- "<="
#'
#' # Now we put together the input P and use it to call GloptiPolyR
#' P <- list(p_f, p_g_1, p_g_2, p_g_3, p_g_4, p_g_5, p_g_6)
#' GloptiPolyR(P)
#' @export
GloptiPolyR <- function(P) {
  # Check polynomial order --------------------------------------------------
  dimensions <- dim(P[[1]]$c)
  degree <- dimensions[1] - 1
  if (degree < 2 || degree > 3) {
    stop("This function accepts only quadratic or cubic polynomials!")
  }
  # Check polynomiao order end ----------------------------------------------
  # Check number of variables -----------------------------------------------
  k <- length(dim(P[[1]]$c))
  if (k < 2 || k > 5) {
    stop("This function accepts only 2 - 5 variables!")
  }
  # Check number of variables end -------------------------------------------
  order <- NULL
  tol <- 1e-4
  maxiter <- 20
  ########################
  ## Secondary functions
  ########################
  # The function 'generate' generates all numbers in base 'base' with 'digits' digits that sum to 'sum', and it returns their
  # equivalent decimal (base 10) representations. For example, generate(3,2,5) produces the vector [2,6,26,10,30,50], which in base 5
  # is [002,011,101,020,110,200].
  generate <- function(digits, sum, base) {
    if (digits < 1) {
      t <- c()
    }
    else if (digits < 2) {
      t <- sum
    }
    else {
      t <- rep(0, choose(sum + digits - 1, digits - 1))
      j <- 1
      for (i in sum:0) {
        s <- generate(digits - 1, sum - i, base)
        r <- length(s)
        t[j:(j + r - 1)] <- i * rep(1, r) + s * base
        j <- j + r
      }
    }
    return(t)
  }
  # The function 'powers2base' converts a vector of exponents of some monomial to a number in base 'base', and it returns its
  # equivalent decimal (base 10) representation. For example, powers2base(c(3,2),5) returns 13. The input vector c(3,2) is meant to
  # represent the exponents of the monomial x[1]^3*x[2]^2.
  powers2base <- function(powers_seq, base) {
    num <- 0
    for (i in 1:length(powers_seq)) {
      num <- num + powers_seq[i] * base^(i - 1)
    }
    return(num)
  }
  # The function 'generate_moment_matrix' generates the moment matrix M_order(y) corresponding to 'nb_var' variables and base 'base'.
  # The entries are the indices of the moments in decimal (base 10) equivalent form. For example, generate_moment_matrix(2,2,5) produces
  #   matrix( c(0,  1,  5,  2,  6,  10),
  #           c(1,  2,  6,  3,  7,  11),
  #           c(5,  6,  10, 7,  11, 15),
  #           c(2,  3,  7,  4,  8,  12),
  #           c(6,  7,  11, 8,  12, 16),
  #           c(10, 11, 15, 12, 16, 20), 6, 6)
  # which in base 5 is
  #   matrix( c(00, 10, 01, 20, 11, 02),
  #           c(10, 20, 11, 30, 21, 12),
  #           c(01, 11, 02, 21, 12, 03),
  #           c(20, 30, 21, 40, 31, 22),
  #           c(11, 21, 12, 31, 22, 13),
  #           c(02, 12, 03, 22, 13, 04), 6, 6).
  generate_moment_matrix <- function(nb_var, order, base) {
    moment_firstrow <- c(0)
    moment_matrix <- matrix(0, 1, 1)
    if (order > 0) {
      for (i in 1:order) {
        moment_firstrow <- c(moment_firstrow, generate(nb_var, i, base))
      }
      moment_matrix <- matrix(0, length(moment_firstrow), length(moment_firstrow))
      for (i in 1:length(moment_firstrow)) {
        moment_matrix[i, ] <- moment_firstrow + moment_firstrow[i]
      }
    }
    return(moment_matrix)
  }
  # The function 'evaluate' evaluates a polynomial, whose coefficients are given in array 'A', at the vector 'x'.
  evaluate <- function(A, x) {
    nonzeros_ind <- which(A != 0)
    powers <- which(A != 0, arr.ind = T)
    powers <- powers - 1
    value <- 0
    for (i in 1:nrow(powers)) {
      value_temp <- 1
      for (j in 1:ncol(powers)) {
        value_temp <- value_temp * x[j]^powers[i, j]
      }
      value <- value + A[nonzeros_ind[i]] * value_temp
    }
    return(as.numeric(value))
  }
  ##################################################
  ## Convert argument P into arguments for Rdsdp
  ##################################################
  # Begin loop that iterates through increasing orders of relaxation until the global
  # ptimum is found (max iterations = 5)
  iter <- 1
  stop <- F
  # tic <- Sys.time()
  while (stop == F & iter < maxiter) {
    # Number of constraints, variables, and basis elements; max degree
    nb_constr <- length(P) - 1 # number of constraints
    nb_var <- 0 # number of variables
    deg <- 0 # max degree over objective and constraints
    for (i in 1:(nb_constr + 1)) {
      # vector
      if (is.null(dim(P[[i]]$c))) {
        nb_var <- 1
        P[[i]]$d <- length(P[[i]]$c) - 1
      }
      # matrix or array
      else {
        nb_var <- max(nb_var, length(dim(P[[i]]$c)))
        powers <- which(P[[i]]$c != 0, arr.ind = T)
        powers <- powers - 1
        P[[i]]$d <- max(rowSums(powers))
      }
      deg <- max(deg, P[[i]]$d)
    }
    # obtain indices of objective function and equality and inequality constraints
    obj_ind <- c()
    eq_ind <- c()
    ineq_ind <- c()
    for (i in 1:(nb_constr + 1)) {
      if (P[[i]]$t == "min" || P[[i]]$t == "max") {
        obj_ind <- c(obj_ind, i)
      }
      else if (P[[i]]$t == "==") {
        eq_ind <- c(eq_ind, i)
      }
      else if (P[[i]]$t == ">=" || P[[i]]$t == "<=") {
        ineq_ind <- c(ineq_ind, i)
      }
    }
    # order of relaxation, base, order of moment matrices
    if (is.null(order)) {
      order <- ceiling(deg / 2) # order of relaxaion; order of main moment matrix
    }
    base <- 2 * order + 1
    for (i in 1:nb_constr + 1) {
      P[[i]]$o <- order - ceiling(P[[i]]$d / 2) # order of moment matrices in LMIs
    }
    P[[obj_ind]]$o <- order # order of moment matrix stored in objective function list
    # change max to min
    if (P[[obj_ind]]$t == "max") {
      P[[obj_ind]]$c <- -P[[obj_ind]]$c
    }
    # change inequality direction to '>='
    for (i in ineq_ind) {
      if (P[[i]]$t == "<=") {
        P[[i]]$c <- -P[[i]]$c
        P[[i]]$t <- ">="
      }
    }
    # sort P so that it has the order: equality constraints, objective function/moment matrix, inequality constraints
    indices <- c(eq_ind, obj_ind, ineq_ind)
    P <- P[indices]
    # new indices
    if (!is.null(eq_ind)) {
      eq_ind <- c(1:length(eq_ind))
    }
    obj_ind <- length(eq_ind) + 1
    ineq_ind <- c((2 + length(eq_ind)):length(P))
    # generate indices for each basis element (excluding the constant term) in base 'base'; var_ind[i] gives the variable with
    # basis index i
    moment_matrix_ind <- generate_moment_matrix(nb_var, order, base)
    basis_ind <- unique(as.vector(moment_matrix_ind))
    basis_ind <- basis_ind[2:length(basis_ind)] # remove 0, the index for the basis element 1
    var_ind <- rep(0, max(basis_ind))
    for (i in 1:length(basis_ind)) {
      var_ind[basis_ind[i]] <- i
    }
    # shift in indices for A and C
    shift <- c(0)
    dims <- c()
    if (!is.null(eq_ind)) {
      for (i in 1:length(eq_ind)) {
        dim_temp <- choose(nb_var + P[[i]]$o, nb_var)
        dim_temp <- (dim_temp + 1) * dim_temp / 2
        dims <- c(dims, dim_temp)
        shift[i + 1] <- dims[i]
      }
    }
    shift[obj_ind + 1] <- choose(nb_var + order, nb_var)^2
    for (i in ineq_ind) {
      shift[i + 1] <- choose(nb_var + P[[i]]$o, nb_var)^2
    }
    shift <- cumsum(shift)
    # Prep b
    b <- rep(0, length(basis_ind))
    nonzeros <- which(P[[obj_ind]]$c != 0)
    powers <- arrayInd(nonzeros, dim(P[[obj_ind]]$c)) - 1
    constant <- 0 # extract constant of objective function
    for (i in 1:nrow(powers)) {
      basis_index <- powers2base(powers[i, ], base)
      var <- var_ind[basis_index]
      if (length(var) == 0) {
        constant <- -P[[obj_ind]]$c[nonzeros[i]]
      }
      b[var] <- -P[[obj_ind]]$c[nonzeros[i]]
    }
    # Prep A & C
    A <- matrix(0, length(basis_ind), shift[length(shift)])
    C <- rep(0, shift[length(shift)])
    # equality constraint
    if (!is.null(eq_ind)) {
      for (i in eq_ind) {
        base_moment_matrix_ind <- generate_moment_matrix(nb_var, P[[i]]$o, base)
        nonzeros <- which(P[[i]]$c != 0)
        powers <- arrayInd(nonzeros, dim(P[[i]]$c)) - 1
        terms <- c()
        for (j in 1:nrow(powers)) {
          terms <- c(terms, powers2base(powers[j, ], base))
        }
        upper_tri_ind <- upper.tri(base_moment_matrix_ind, diag = T)
        upper_tri_ind <- which(upper_tri_ind == T)
        for (j in 1:dims[i]) {
          for (k in 1:length(terms)) {
            var <- var_ind[base_moment_matrix_ind[upper_tri_ind[j]] + terms[k]]
            if (length(var) == 0) {
              C[shift[i] + j] <- P[[i]]$c[nonzeros[k]]
            }
            else {
              A[var, shift[i] + j] <- -P[[i]]$c[nonzeros[k]]
            }
          }
        }
      }
    }
    # moment matrix
    for (i in 1:length(moment_matrix_ind)) {
      var <- var_ind[moment_matrix_ind[i]]
      if (moment_matrix_ind[i] == 0) {
        C[shift[obj_ind] + i] <- 1
      }
      else {
        A[var, shift[obj_ind] + i] <- -1
      }
    }
    dims <- c(dims, nrow(moment_matrix_ind)) # dimensions of LMIs
    # LMIs for each inequality constraint
    for (i in ineq_ind) {
      base_moment_matrix_ind <- generate_moment_matrix(nb_var, P[[i]]$o, base)
      dims <- c(dims, nrow(base_moment_matrix_ind))
      nonzeros <- which(P[[i]]$c != 0)
      powers <- arrayInd(nonzeros, dim(P[[i]]$c)) - 1
      terms <- c()
      for (j in 1:nrow(powers)) {
        terms <- c(terms, powers2base(powers[j, ], base))
      }
      for (j in 1:length(base_moment_matrix_ind)) {
        for (k in 1:length(terms)) {
          var <- var_ind[base_moment_matrix_ind[j] + terms[k]]
          if (length(var) == 0) {
            C[shift[i] + j] <- P[[i]]$c[nonzeros[k]]
          } else {
            A[var, shift[i] + j] <- -P[[i]]$c[nonzeros[k]]
          }
        }
      }
    }
    # Prep K
    K <- NULL
    K$l <- sum(dims[eq_ind])
    K$s <- c(dims[obj_ind], dims[ineq_ind])
    ################
    ## Call Rdsdp
    ################
    # Call Rdsdp to solve LMI relaxation
    OPTIONS <- list()
    OPTIONS$print <- 0
    results <- Rdsdp::dsdp(A, b, C, K, OPTIONS)
    # optimal value of objective function for relaxed problem
    Qr_obj <- results$STATS$dobj + constant
    if (P[[obj_ind]]$t == "min") {
      Qr_obj <- -Qr_obj
    }
    #########################################
    ## Check if optimality has been reached
    #########################################
    # Check if LMI relaxed vector is the optimal solution to (P)
    x <- results$y[1:nb_var]
    value <- evaluate(P[[obj_ind]]$c, x)
    dist <- abs(value - Qr_obj)
    if (dist < tol) {
      feasible <- T
      i <- 1
      while (feasible & i <= (nb_constr + 1)) {
        value <- evaluate(P[[i]]$c, x)
        if (P[[i]]$t == ">=" & value < (-tol)) {
          feasible <- F
        } else if (P[[i]]$t == "==" & abs(value) > tol) {
          feasible <- F
        }
        i <- i + 1
      }
      if (feasible) {
        stop <- T
        # cat("Global optimum has been found.", "\n")
        # cat("Order of relaxation: ", order, "\n")
        # cat("x* = ", round(x,4), "\n")
        # cat("p* = ", round(Qr_obj,4), "\n")
      } else {
        # cat("Global optimum has not been found. Increasing the order of relaxation from ", order, "to ", order+1, ".", "\n")
      }
    } else {
      # cat("Global optimum has not been found. Increasing the order of relaxation from ", order, "to ", order+1, ".", "\n")
    }
    # update iteration counter; end of while loop
    iter <- iter + 1
    order <- order + 1
    if (order > 5) stop("gloptipoly takes too many relaxation steps")
    # when set order = 10, n = 200, memory exhausted error
    # toc <- Sys.time()
    # if((toc - tic) > 5) stop("gloptipoly solver takes too long")
  }
  # Extract moment matrices of increasing order up to 'order'
  # psd_constr = C-t(A)%*%results$y
  # moment_matrix = matrix(psd_constr[(shift[obj_ind]+1):(shift[obj_ind]+dims[obj_ind]^2)],dims[obj_ind],dims[obj_ind])
  # moment_matrices = list()
  # for(i in 1:order){
  #  dim_temp = choose(nb_var+i,nb_var)
  #  moment_matrices[[i]] = moment_matrix[1:dim_temp,1:dim_temp]
  # }
  ###############
  ## Extract x*
  ###############
  # end gloptipoly
  # return
  list(solution = round(x, 4), objective = round(Qr_obj, 4))
}
