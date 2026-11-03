#' HDLCP: High-dimensional mean change-point inference
#' @keywords internal
#' @import stats
#' @importFrom utils write.csv read.table download.file unzip head tail
#' @import graphics
#' @importFrom grDevices pdf dev.off
#' @useDynLib HDLCP, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

.integer <- function(x, name, lower = 1L) {
  if (length(x) != 1L || !is.finite(x) || x != floor(x) || x < lower)
    stop(name, " must be an integer >= ", lower, call. = FALSE)
  as.integer(x)
}

.data <- function(X, check_variance=TRUE) {
  if (!is.matrix(X) || !is.numeric(X) || nrow(X) < 25L || ncol(X) < 2L ||
      any(!is.finite(X))) stop("X must be a finite numeric matrix with n >= 25 and p >= 2.")
  storage.mode(X) <- "double"
  if (check_variance && any(colSums((X[-1, , drop = FALSE] - X[-nrow(X), , drop = FALSE])^2) <= 0))
    stop("Every column must have positive adjacent-difference variance.")
  X
}

.with_seed <- function(seed, fun) {
  old_kind <- RNGkind()
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", .GlobalEnv)
  on.exit({do.call(RNGkind, as.list(old_kind)); if (had) assign(".Random.seed", old, .GlobalEnv)
           else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)})
  set.seed(.integer(seed, "seed", 0L))
  fun()
}

.parallel <- function(tasks, fun, workers, seed) {
  workers <- min(.integer(workers, "workers"), length(tasks))
  seeds <- .with_seed(seed, function() sample.int(.Machine$integer.max, length(tasks)))
  job <- function(i) .with_seed(seeds[i], function() fun(tasks[[i]]))
  if (workers == 1L) return(lapply(seq_along(tasks), job))
  cl <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterCall(cl, function(paths) {
    .libPaths(paths)
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    loadNamespace("HDLCP")
    NULL
  }, .libPaths())
  parallel::parLapply(cl, seq_along(tasks), job)
}

.blocks <- function(m, count) {
  sizes <- rep((m - count + 1L) %/% count, count)
  extra <- (m - count + 1L) %% count
  if (extra) sizes[seq_len(extra)] <- sizes[seq_len(extra)] + 1L
  start <- cumsum(c(1L, head(sizes + 1L, -1L)))
  Map(function(a, b) seq.int(a, length.out = b), start, sizes)
}

.prepare <- function(X, lambda) {
  n <- nrow(X); p <- ncol(X); r <- seq_len(n - 1L); t <- r/n; v <- t*(1-t)
  D <- (X[-1, , drop=FALSE] - X[-n, , drop=FALSE])/sqrt(2)
  sigma2 <- colMeans(D^2)
  cs <- apply(X, 2, cumsum)
  W <- sweep((cs[-n, , drop=FALSE] - tcrossprod(t, cs[n, ]))/sqrt(n), 2, sqrt(sigma2), "/")
  ordered <- t(apply(W^2, 1, sort, decreasing=TRUE))
  list(n=n, p=p, v=v, trim=which(r >= lambda & r <= n-lambda),
       Tn=2*log((n-lambda)/lambda), sigma2=sigma2,
       D=sweep(D, 2, sqrt(sigma2), "/"),
       sums=t(apply(ordered, 1, cumsum)), W=W)
}

.summaries <- function(prepared) {
  D <- prepared$D; m <- nrow(D); p <- ncol(D)
    LR <- max(1L, min(ceiling(16*log((m+1)*p)), floor((m+1)/2)))
    blocks <- .blocks(m, LR)
    entries <- vapply(blocks, function(i) as.vector(crossprod(D[i,,drop=FALSE])/length(i)), numeric(p*p))
    R <- matrix(apply(entries, 1, median), p, p)
    traces <- vapply(.blocks(m, 5L), function(i) {
      G <- tcrossprod(D[i,,drop=FALSE]); use <- abs(outer(i,i,"-")) > 1
      mean(G[use]^2)
    }, numeric(1))
    trace <- median(traces)
  diag(R) <- 1
  list(R=R, s2=max(trace/p-1,0), threshold=2*sqrt(log(p)/m))
}

#' Combine component p-values using the Cauchy rule
#' @param p Numeric vector of p-values in [0,1].
#' @param weights Positive weights, normalized to sum to one.
#' @return A combined p-value, using the endpoint extension in Appendix A:
#' any zero yields zero; otherwise any one yields one.
#' @export
#' @examples
#' cauchy_combine(c(0.01, 0.2, 0.5))
cauchy_combine <- function(p, weights = rep(1/length(p), length(p))) {
  if (!length(p) || any(!is.finite(p)) || any(p<0|p>1) || length(weights)!=length(p) ||
      any(!is.finite(weights)|weights<=0)) stop("Invalid p-values or weights.")
  if (any(p==0)) return(0)
  if (any(p==1)) return(1)
  w <- weights/sum(weights)
  terms <- ifelse(p < 0.5, 1/tan(pi*p), -1/tan(pi*(1-p)))
  z <- sum(w*terms)
  if (z > 0) atan(1/z)/pi else 1-atan(-1/z)/pi
}
