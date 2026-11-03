#' Test and locate a single high-dimensional mean change
#' @param X Numeric n by p matrix; rows ordered observations, columns variables.
#' @param reference Reusable object from [hdlcp_calibrate()].
#' @param alpha Test level in (0,1).
#' @param gamma Boundary weights, any subset of c(0,0.5).
#' @param keep_paths Whether to return component scan paths.
#' @return An hdlcp_fit containing component and Cauchy p-values, rejections,
#' location estimates, and data-driven nuisance summaries. Locations are the
#' last row before the change. Cauchy localization selects the component with
#' the smallest global p-value; exact ties use component order (fixed, then proportional).
#' The log.p.value column preserves analytic tail comparisons under numerical
#' underflow. Reference-bank endpoint values are flagged in diagnostics.
#' @details The reference must match the dimensions of X. Neither the true
#' covariance, change point, nor innovation law is accepted by this function.
#' Four-component Cauchy calibration is assessed empirically in v37; the
#' paper's independence theorem establishes a matched fixed/proportional pair.
#' @export
#' @examples
#' \donttest{
#' ref <- hdlcp_calibrate(40,20,k=1,q=.25,B=99)
#' X <- simulate_hdlcp_data(40,20,s=3,delta=1)
#' fit <- hdlcp_test(X,ref)
#' fit$tests
#' }
hdlcp_test <- function(X, reference, alpha=.05, gamma=c(0,.5), keep_paths=FALSE) {
  X <- .data(X,check_variance=FALSE); ref <- reference
  if(!inherits(ref,"hdlcp_reference") || nrow(X)!=ref$n || ncol(X)!=ref$p)
    stop("Reference and data dimensions do not match.")
  if(!is.null(ref$version) && utils::compareVersion(ref$version,"0.1.2")<0)
    stop("Rebuild this numerical reference with HDLCP >= 0.1.2; its tangent tables are incompatible.")
  if(length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=1)stop("Invalid alpha.")
  if(!length(gamma)||any(!gamma %in% c(0,.5))||anyDuplicated(gamma))stop("Invalid gamma.")
  if(any(colSums(diff(X)^2)==0)) {
    rows <- do.call(rbind,lapply(gamma,function(g) {
      names <- c(paste0("F",ref$k,"-",g)[seq_along(ref$k)],
                 paste0("P",ref$q,"-",g)[seq_along(ref$q)],paste0("CL",g))
      data.frame(method=names,gamma=g,statistic=0,p.value=1,log.p.value=0,reject=FALSE,
        location=if(g==0)1L else ref$lambda,selected=names)
    }))
    return(structure(list(tests=rows,components=rows[!grepl("^CL",rows$method),],
      diagnostics=list(degenerate_scale=TRUE),paths=list(),alpha=alpha,n=ref$n,p=ref$p,
      reference=list(B=ref$B,seed=ref$seed)),class="hdlcp_fit"))
  }
  pr <- .prepare(X,ref$lambda); sm <- .summaries(pr)
  rows <- list(); paths <- list(); diag <- list(s2=sm$s2,sigma2=pr$sigma2)
  add <- function(method,g,stat,pv,loc,path,reject=pv<alpha,logpv=log(pv)) {
    rows[[length(rows)+1L]] <<- data.frame(method=method,gamma=g,statistic=stat,
      p.value=pv,log.p.value=logpv,reject=reject,location=as.integer(loc),selected=method)
    if(g==0 || startsWith(method,"F"))
      diag$reference_boundary[[method]] <<- pv==0 || pv==1
    if(keep_paths) {
      if(g==.5)path[-pr$trim]<-NA_real_
      paths[[method]] <<- path
    }
  }
  for(k in ref$k) {
    th <- .theta(ref,k,sm); diag[[paste0("theta_k",k)]] <- th
    for(g in gamma) {
      path <- pr$sums[,k]/pr$v^(2*g); allowed <- if(g==0)seq_along(path) else pr$trim
      loc <- allowed[which.max(path[allowed])]; stat <- path[loc]
      cv <- .fixed_curve(ref,k,th,g)
      pv <- .curve_p(stat,cv,ref$tails)
      add(paste0("F",k,"-",g),g,stat,pv,loc,path,
          reject=stat>.critical_at(alpha,cv,ref$tails))
    }
  }
  m <- ref$n-1; kap <- (3*m-1)/m^2
  ell <- log(ref$n*ref$p); rate <- ref$p^(-.5)+ref$n^(-.5)+sqrt(ell/ref$n)+ell^2/ref$n
  for(j in seq_along(ref$q)) {
    q <- ref$q[j]; lu <- ref$lookups[[j]]; kp <- ceiling(q*ref$p)
    ker <- .kernel(pr,sm,lu,need_factor=0 %in% gamma)
    center <- ref$p*lu$mu+(kp-ref$p*q)*lu$uq+ref$p*kap*(lu$mu+.5*lu$uq^2*dchisq(lu$uq,1))
    U <- rowSums(pmax(pr$D^2-lu$uq,0)+q*lu$uq-lu$mu)/sqrt(ref$p)
    V <- U-mean(U)
    raw <- sqrt(m*(m-1))/(m-2)*mean(V^3)/(mean(V^2)+rate^2)^1.5
    skew <- max(-rate,min(rate,raw))
    diag[[paste0("q",q)]] <- c(omega2=ker$omega,eta=ker$eta,skewness=skew,strong=ker$strong)
    for(g in gamma) {
      path <- (pr$sums[,kp]-pr$v*center)/(sqrt(ref$p*ker$omega)*pr$v^(2*g))
      allowed <- if(g==0)seq_along(path) else pr$trim
      loc <- allowed[which.max(path[allowed])]; stat <- path[loc]
      if(g==0) {
        G <- ker$factor %*% ref$normal
        maxima <- apply(pr$v*(G+skew*(G^2-1)/6),2,max)
        pv <- .reference_tail(stat,maxima)
        logpv <- log(pv)
      } else {
        z <- sqrt(2*log(pr$Tn))*stat - (2*log(pr$Tn)+.5*log(log(pr$Tn))+log(ker$eta/sqrt(pi)))
        logpv <- .gumbel_log_tail(z)
        pv <- exp(logpv)
      }
      add(paste0("P",q,"-",g),g,stat,pv,loc,path,logpv=logpv)
    }
  }
  components <- do.call(rbind,rows)
  combined <- lapply(gamma,function(g) {
    z <- components[components$gamma==g,,drop=FALSE]; best <- which.min(z$log.p.value)
    cc <- .cauchy_result(z$p.value)
    data.frame(method=paste0("CL",g),gamma=g,statistic=cc$statistic,p.value=cc$p.value,
      log.p.value=log(cc$p.value),
      reject=cc$p.value<alpha,location=z$location[best],selected=z$method[best])
  })
  structure(list(tests=rbind(components,do.call(rbind,combined)),components=components,
    diagnostics=diag,paths=paths,alpha=alpha,n=ref$n,p=ref$p,
    reference=list(B=ref$B,seed=ref$seed,floor_level=ref$floor_level,
      tail_resolution=1/ref$B)),class="hdlcp_fit")
}

#' @export
print.hdlcp_fit <- function(x, ...) {
  cat("HDLCP single-change inference: n =",x$n,", p =",x$p,"\n")
  print(x$tests,row.names=FALSE,...)
  invisible(x)
}

#' Extract component-selected change-point locations
#' @param object An hdlcp_fit.
#' @return Named vector of location estimates (last row before the change).
#' @export
hdlcp_location <- function(object) {
  if(!inherits(object,"hdlcp_fit"))stop("Expected an hdlcp_fit.")
  setNames(object$tests$location,object$tests$method)
}

#' Evaluate one fixed-order component
#' @param X Numeric observation matrix.
#' @param reference Reusable reference containing k.
#' @param k Fixed order.
#' @param gamma Boundary weight, zero or one half.
#' @param alpha Test level.
#' @return One-row data frame with statistic, p.value, rejection and location.
#' @export
hdlcp_fixed_test <- function(X,reference,k=1L,gamma=0,alpha=.05) {
  if(!inherits(reference,"hdlcp_reference")||length(k)!=1||!k%in%reference$k)
    stop("k must occur in the reference.")
  reference$k<-k;reference$q<-numeric();reference$lookups<-list()
  hdlcp_test(X,reference,alpha=alpha,gamma=gamma)$components
}

#' Evaluate one proportional-order component
#' @param X Numeric observation matrix.
#' @param reference Reusable reference containing q.
#' @param q Proportion strictly between zero and one.
#' @param gamma Boundary weight, zero or one half.
#' @param alpha Test level.
#' @return One-row data frame with statistic, p.value, rejection and location.
#' @export
hdlcp_proportional_test <- function(X,reference,q=.1,gamma=0,alpha=.05) {
  if(!inherits(reference,"hdlcp_reference")||length(q)!=1||!q%in%reference$q)
    stop("q must occur in the reference.")
  i<-match(q,reference$q);reference$k<-integer();reference$q<-q
  reference$lookups<-reference$lookups[i]
  hdlcp_test(X,reference,alpha=alpha,gamma=gamma)$components
}

#' Plot component scans
#' @param x An hdlcp_fit computed with keep_paths=TRUE.
#' @param ... Arguments passed to matplot.
#' @export
plot.hdlcp_fit <- function(x, ...) {
  if(!length(x$paths))stop("Use keep_paths=TRUE in hdlcp_test().")
  matplot(seq_len(x$n-1),do.call(cbind,x$paths),type="l",lty=1,
          xlab="Candidate split",ylab="Component scan",...)
  legend("topright",names(x$paths),col=seq_along(x$paths),lty=1,cex=.7)
  invisible(x)
}
