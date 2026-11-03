# Numerical evaluation of the Gaussian functions defined in manuscript v37.
gauss_legendre <- function(n, a, b) {
  i <- seq_len(n - 1L)
  off <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- off
  J[cbind(i + 1L, i)] <- off
  ee <- eigen(J, symmetric = TRUE)
  ord <- order(ee$values)
  list(
    x = (a + b) / 2 + (b - a) * ee$values[ord] / 2,
    w = (b - a) * ee$vectors[1L, ord]^2
  )
}

build_q_lookup <- function(q, n_corr = 4001L, n_quad = 80L) {
  if(q == 1) {
    corr <- seq(0,1,length.out=n_corr)
    return(list(q=1,aq=0,uq=0,mu=1,tau2=2,corr=corr,gamma=2*corr^2,gamma_prime=4*corr))
  }
  aq <- qnorm(log(q)-log(2), lower.tail=FALSE, log.p=TRUE)
  uq <- aq^2
  mu <- q + 2 * aq * dnorm(aq)
  dq <- 2 * aq * dnorm(aq) + q * (1 - aq^2)
  eq <- 2 * aq * (3 - aq^2) * dnorm(aq) +
    q * (aq^4 - 2 * aq^2 + 3)
  tau2 <- eq - dq^2

  bound <- max(10, aq+8)
  rules <- list(
    gauss_legendre(n_quad, -bound, -aq),
    gauss_legendre(n_quad, -aq, aq),
    gauss_legendre(n_quad, aq, bound)
  )
  x <- unlist(lapply(rules, `[[`, "x"), use.names = FALSE)
  w <- unlist(lapply(rules, `[[`, "w"), use.names = FALSE) * dnorm(x)
  psi_x <- pmax(x^2 - uq, 0) + q * uq
  quadrature_mean <- sum(w * psi_x)

  corr <- seq(0, 1, length.out = n_corr)
  gamma <- numeric(n_corr)
  gamma_prime <- numeric(n_corr)
  for (ii in seq_len(n_corr - 1L)) {
    a <- corr[ii]
    conditional_sd <- sqrt(pmax(1 - a^2, .Machine$double.eps))
    conditional_mean <- a * x
    upper_z <- (aq - conditional_mean) / conditional_sd
    lower_z <- (-aq - conditional_mean) / conditional_sd
    upper_prob <- pnorm(upper_z, lower.tail = FALSE)
    lower_prob <- pnorm(lower_z)
    tail_prob <- upper_prob + lower_prob
    upper_second <- (conditional_mean^2 + conditional_sd^2) * upper_prob +
      conditional_sd * (conditional_mean + aq) * dnorm(upper_z)
    lower_second <- (conditional_mean^2 + conditional_sd^2) * lower_prob +
      conditional_sd * (aq - conditional_mean) * dnorm(lower_z)
    conditional_psi <- upper_second + lower_second - uq * tail_prob + q * uq
    gamma[ii] <- sum(w * psi_x * conditional_psi) - quadrature_mean * mu
    # Gamma'(rho) = E[h'(Z1) h'(Z2)], h(z)=psi_q(z^2).
    conditional_derivative <- 2*(conditional_mean*tail_prob +
      conditional_sd*(dnorm(upper_z)-dnorm(lower_z)))
    gamma_prime[ii] <- sum(w * (2*x*(abs(x)>aq)) * conditional_derivative)
  }
  gamma[1L] <- 0
  gamma[n_corr] <- tau2
  gamma <- pmax(gamma, 0)
  gamma_prime[1L] <- 0
  gamma_prime[n_corr] <- 4 * mu
  list(
    q = q, aq = aq, uq = uq, mu = mu, tau2 = tau2,
    corr = corr, gamma = gamma, gamma_prime = pmax(gamma_prime, 0)
  )
}

gamma_at <- function(lookup, a) {
  ans <- approx(
    lookup$corr, lookup$gamma, xout = as.vector(abs(a)),
    rule = 2, ties = "ordered"
  )$y
  if (is.null(dim(a))) ans else array(ans, dim = dim(a))
}

gamma_prime_at <- function(lookup, a) {
  ans <- approx(
    lookup$corr, lookup$gamma_prime, xout = as.vector(abs(a)),
    rule = 2, ties = "ordered"
  )$y * sign(as.vector(a))
  if (is.null(dim(a))) ans else array(ans, dim = dim(a))
}

fixed_marginal_tail <- function(x, kappa) {
  x <- pmax(x, 0)
  ans <- pchisq(x, df=1, lower.tail=FALSE)
  use <- is.finite(x) & x>0
  ans[use] <- ans[use]+kappa*x[use]*(x[use]+1)*dchisq(x[use],1)/4
  ans
}

fixed_marginal_density <- function(x, kappa) {
  correction <- 1 + kappa * (x^2 - 2 * x - 1) / 8
  ans <- dchisq(x, df=1)*correction
  ans[is.infinite(x) & x>0] <- 0
  ans
}

.fixed_level <- function(u, kappa) {
  uniroot(function(x) fixed_marginal_tail(x,kappa)-u,
    c(0,max(40,10*log(1/u))),tol=1e-12*min(1,(1-u)^2))$root
}

bivariate_abs_tail <- function(rho, threshold) {
  rho <- min(max(abs(rho), 0), 0.999)
  a <- sqrt(threshold)
  sd_cond <- sqrt(1 - rho^2)
  integrand <- function(z) {
    dnorm(z) * (
      pnorm((rho * z - a) / sd_cond) +
        pnorm((-rho * z - a) / sd_cond)
    )
  }
  2 * integrate(
    integrand, lower = a, upper = Inf,
    rel.tol = 2e-8, subdivisions = 200L
  )$value
}

build_fixed_cluster_lookup <- function(k, p, kappa, rho_grid) {
  lapply(seq_len(k), function(rank) {
    marginal_tail <- rank / p
    threshold <- .fixed_level(marginal_tail,kappa)
    list(
      rank = rank,
      marginal_tail = marginal_tail,
      threshold = threshold,
      pair_tail_grid = vapply(
        rho_grid, bivariate_abs_tail, numeric(1), threshold = threshold
      )
    )
  })
}


