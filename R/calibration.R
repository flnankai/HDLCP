#' Build reusable numerical reference tables
#'
#' All references are independent of the observations and of their covariance.
#' Fixed-order tables combine a finite independent Gaussian anchor and a
#' compound-Poisson tangent increment. Proportional orders reuse standard normal
#' vectors, transformed by the estimated full kernel at test time.
#' No reference is persisted automatically.
#' @param n,p Sample size and number of coordinates.
#' @param k Fixed orders, integer vector strictly smaller than p; may be empty.
#' @param q Proportions strictly between zero and one; may be empty.
#' @param B Number of reference draws, at least 39 (use 2000 or more for inference).
#' @param lambda Boundary trimming, default ceiling(n^0.3).
#' @param theta_grid Cluster grid including zero and 0.30.
#' @param floor_level Tangent particle truncation, negative. Check sensitivity for publication.
#' @param workers Number of parallel workers (Windows supported).
#' @param seed Integer seed; independent of worker count and restores caller RNG.
#' @return An hdlcp_reference object. It can be reused across samples with identical
#' n,p,k,q and lambda. Increasing B controls numerical reference error.
#' @export
#' @examples
#' \donttest{
#' ref <- hdlcp_calibrate(40, 20, k=1, q=0.25, B=99, workers=1)
#' }
hdlcp_calibrate <- function(n, p, k=c(1L,5L), q=c(0.1,0.9), B=2000L,
                            lambda=ceiling(n^0.3),
                            theta_grid=c(0,.025,.05,.075,.10,.15,.20,.30),
                            floor_level=-5, workers=1L, seed=20261103L) {
  n <- .integer(n,"n",25); p <- .integer(p,"p",2); B <- .integer(B,"B",39)
  lambda <- .integer(lambda,"lambda")
  if (lambda>=n/2 || 2*log((n-lambda)/lambda)<=1) stop("lambda gives an inadmissible logit interval.")
  if (any(!is.finite(k)|k!=floor(k)|k<1|k>=p) || anyDuplicated(k)) stop("Invalid k.")
  if (any(!is.finite(q)|q<=0|q>=1) || anyDuplicated(q)) stop("q must be distinct and in (0,1).")
  if (!length(k)&&!length(q)) stop("Provide k or q.")
  if (any(!is.finite(theta_grid)) || length(theta_grid)<2 || min(theta_grid)!=0 ||
      max(theta_grid)!=.3 || any(diff(theta_grid)<=0)) stop("Invalid theta grid.")
  if (length(floor_level)!=1 || !is.finite(floor_level) || floor_level>=0) stop("floor_level must be negative.")
  m <- n-1L; kap <- (3*m-1)/m^2; t <- seq_len(m)/n; v <- t*(1-t)
  trim <- which(seq_len(m)>=lambda & seq_len(m)<=n-lambda)
  tails <- sort(unique(c(seq_len(B)/(B+1),.01,.025,.05,.1)))
  tails <- tails[tails>=1/(B+1) & tails<=B/(B+1)]
  curves <- function(z) vapply(seq_len(ncol(z)),function(j)
    as.numeric(quantile(z[,j],1-tails,type=8)),numeric(length(tails)))
  fixed <- list()
  if (length(k)) {
    groups <- split(seq_len(B), ceiling(seq_len(B)/32))
    anchor <- do.call(rbind,.parallel(groups,function(ids) {
      t(vapply(ids,function(i) {
        pr <- .prepare(matrix(rnorm(n*p),n,p),lambda)
        unlist(lapply(k,function(kk)c(max(pr$sums[,kk]),max(pr$sums[trim,kk]/v[trim]))))
      },numeric(2*length(k))))
    },workers,seed))
    beta <- uniroot(function(x) fixed_marginal_tail(x,kap)-1/p,c(1e-8,max(40,10*log(p))))$root
    scale <- fixed_marginal_tail(beta,kap)/fixed_marginal_density(beta,kap)
    tasks <- expand.grid(ki=seq_along(k),ti=seq_along(theta_grid))
    banks <- .parallel(split(tasks,seq_len(nrow(tasks))),function(task) {
      z <- tangent_bank_cpp(B,theta_grid[task$ti],k[task$ki],beta,scale,
                            beta*qlogis(t),v,trim,floor_level)
      if(any(!is.finite(z))) stop("Tangent truncation too shallow; decrease floor_level.")
      curves(z)
    },workers,seed+1L)
    rho_grid <- seq(0,.95,.005)
    for (i in seq_along(k)) {
      indices <- which(tasks$ki==i)
      tang <- array(unlist(banks[indices]),c(length(tails),2,length(theta_grid)))
      fixed[[as.character(k[i])]] <- list(anchor=curves(anchor[,c(2*i-1,2*i),drop=FALSE]),
        tangent=tang, lookup=build_fixed_cluster_lookup(k[i],p,kap,rho_grid))
    }
  }
  noise <- {
    LR <- max(1,min(ceiling(16*log(n*p)),floor((m+1)/2)))
    .with_seed(seed+2L,function() mean(replicate(B, {
      Y <- matrix(rnorm(2*n),n,2); D <- diff(Y)/sqrt(2)
      D <- sweep(D,2,sqrt(colMeans(D^2)),"/")
      median(vapply(.blocks(m,LR),function(i)mean(D[i,1]*D[i,2]),numeric(1)))^2
    })))
  }
  lookups <- lapply(q, build_q_lookup)
  normal <- if(length(q)) .with_seed(seed+3L,function()matrix(rnorm(m*B),m,B)) else NULL
  structure(list(n=n,p=p,k=k,q=q,B=B,lambda=lambda,theta_grid=theta_grid,
    floor_level=floor_level,seed=seed,tails=tails,fixed=fixed,
    lookups=lookups,normal=normal,noise=noise,version="0.1.1"),class="hdlcp_reference")
}

.fixed_curve <- function(ref,k,theta,gamma) {
  x <- ref$fixed[[as.character(k)]]; j <- if(gamma==0)1L else 2L
  i <- min(findInterval(theta,ref$theta_grid),length(ref$theta_grid)-1L); i <- max(i,1L)
  w <- (theta-ref$theta_grid[i])/diff(ref$theta_grid[c(i,i+1L)])
  raw <- x$anchor[,j]+(1-w)*x$tangent[,j,i]+w*x$tangent[,j,i+1]-x$tangent[,j,1]
  rev(cummax(rev(raw)))
}

.curve_p <- function(stat,crit,tails) {
  # Generalized inverse on a linearly interpolated decreasing curve.
  if(stat<=min(crit))return(tail(tails,1))
  if(stat>max(crit))return(tails[1])
  i <- which(stat>crit)[1]
  if(i==1)return(tails[1])
  if(crit[i-1]==crit[i])return(tails[i])
  tails[i-1]+(crit[i-1]-stat)/(crit[i-1]-crit[i])*(tails[i]-tails[i-1])
}

.theta <- function(ref,k,summ) {
  rho <- abs(summ$R[upper.tri(summ$R)]); strong <- rho[rho>summ$threshold]
  weak <- max(ref$p*summ$s2/2-sum(pmax(strong^2-ref$noise,0)),0)
  mean(vapply(ref$fixed[[as.character(k)]]$lookup,function(lu) {
    u <- lu$marginal_tail; ug <- pchisq(lu$threshold,1,lower.tail=FALSE)
    inflation <- (u/ug)^2
    pair <- approx(seq(0,.95,.005),lu$pair_tail_grid,pmin(strong,.95),rule=2)$y
    es <- sum(pmax(pmin(pair*inflation,u)-u^2,0))
    ew <- 2*lu$threshold*dnorm(sqrt(lu$threshold))^2*inflation*weak
    A <- 2*(es+ew)/(ref$p*u)
    min(.3,A/(2+A))
  },numeric(1)))
}

.kernel <- function(pr,summ,lookup) {
  p <- pr$p; t <- seq_len(pr$n-1)/pr$n; v <- pr$v
  a <- (outer(t,t,pmin)-tcrossprod(t))/sqrt(outer(v,v))
  rho <- pmin(abs(summ$R[upper.tri(summ$R)]),1)
  rho <- rho[rho>summ$threshold]
  C <- gamma_at(lookup,a)+2*lookup$mu^2*summ$s2*a^2
  omega <- lookup$tau2+2*lookup$mu^2*summ$s2
  slope <- 2*lookup$mu+2*lookup$mu^2*summ$s2
  # Full Gamma and Gamma-prime corrections for every selected strong pair.
  if(length(rho)) {
    for(r in rho) C <- C+2/p*(gamma_at(lookup,r*a)-2*lookup$mu^2*r^2*a^2)
    omega <- omega+2/p*sum(gamma_at(lookup,rho)-2*lookup$mu^2*rho^2)
    slope <- slope+1/p*sum(rho*gamma_prime_at(lookup,rho)-4*lookup$mu^2*rho^2)
  }
  if(omega<=0 || slope<=0)stop("Nonpositive estimated variance or intensity.")
  K <- (C+t(C))/(2*omega); diag(K)<-1
  ee <- eigen(K,symmetric=TRUE)
  if(min(ee$values) < -1e-7) {
    if(!requireNamespace("CVXR",quietly=TRUE))
      stop("Kernel needs the v37 max-norm correlation projection; install CVXR.")
    K <- .project_kernel(K); ee <- eigen((K+t(K))/2,symmetric=TRUE)
  }
  F <- sweep(ee$vectors,2,sqrt(pmax(ee$values,0)),"*")
  list(factor=F,omega=omega,eta=slope/omega,strong=length(rho))
}

.project_kernel <- function(K) {
  A <- CVXR::Variable(c(nrow(K),ncol(K)),symmetric=TRUE)
  problem <- CVXR::Problem(CVXR::Minimize(CVXR::max_entries(abs(A-K))),
    constraints=list(CVXR::diag(A)==1,CVXR::`%>>%`(A,0)))
  CVXR::psolve(problem,feastol=1e-9,reltol=1e-9,abstol=1e-9)
  if(!CVXR::status(problem) %in% c("optimal","optimal_inaccurate"))
    stop("Max-norm correlation projection failed.")
  projected <- CVXR::value(A)
  if(any(!is.finite(projected)) || min(eigen(projected,symmetric=TRUE,only.values=TRUE)$values) < -1e-6)
    stop("Correlation projection returned an invalid matrix.")
  projected
}
