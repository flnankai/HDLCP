.dms_required <- function() {
  if(!requireNamespace("DMS",quietly=TRUE))
    stop("Install the published comparator package with remotes::install_github('ghwang-nk/DMS').")
}

#' Wang-Feng DMS test and minimum-p component localization
#' @param X Numeric n by p matrix.
#' @param gamma Boundary weighting, zero or one half.
#' @param lambda Trim for the weighted maximum, default floor(0.1*n).
#' @param type Combination rule accepted by DMS, default Fisher.
#' @return List with p.value, location, and component p-values.
#' @export
dms_test <- function(X,gamma=0,lambda=floor(.1*nrow(X)),type="Fisher") {
  X <- .data(X); .dms_required()
  if(!gamma %in% c(0,.5)||length(gamma)!=1)stop("Invalid gamma.")
  lambda <- .integer(lambda,"lambda")
  if(lambda>=nrow(X)/2)stop("lambda must be below n/2.")
  fit <- DMS::DMS(X,gam=gamma,lam=lambda,type=type)
  paths <- DMS::cusum(X,c(0,.5)); se <- sqrt(drop(DMS::est_var_diff(X)$s2))
  gindex <- if(gamma==0)1L else 2L
  M <- sweep(paths[,,gindex],2,se,"/"); allowed <- if(gamma==0)seq_len(nrow(X)-1) else lambda:(nrow(X)-lambda)
  locmax <- allowed[which.max(apply(abs(M[allowed,,drop=FALSE]),1,max))]
  sumfit <- DMS::WZWY19(X,CUSUM=paths[,,2],PE=FALSE)
  locsum <- which.max(rowSums(sumfit$CUSUM2_s))
  list(p.value=fit$pv.DMS,location=if(fit$pv.Max<=fit$pv.Sum)locmax else locsum,
       component.p=c(max=fit$pv.Max,sum=fit$pv.Sum))
}

#' Zhang-Wang-Shao self-normalized adaptive test
#' @param X Numeric n by p matrix.
#' @return List with p.value, selected location, and q=2,6 component results.
#' @details Uses the algebraically accelerated project implementation and the
#' original DMS calibration samples. The location maximizes the self-normalized
#' scan of the component with the smallest p-value.
#' @export
zws22_test <- function(X) {
  X <- .data(X); .dms_required()
  z <- zws22_fast_stats(X)
  ns <- asNamespace("DMS")
  pv <- c(mean(get("cal_q2",ns)>=z$ts2),mean(get("cal_q6",ns)>=z$ts6))
  loc <- c(3L+which.max(z$num2/pmax(z$denom2,1e-300)),
           11L+which.max(z$num6/pmax(z$denom6,1e-300)))
  list(p.value=1-(1-min(pv))^2,location=loc[which.min(pv)],
       component.p=setNames(pv,c("q2","q6")),statistic=c(z$ts2,z$ts6))
}

#' Liu-Zhou-Zhang-Liu adaptive norm test
#' @param X Numeric n by p matrix.
#' @param B Number of Gaussian multiplier draws.
#' @param boundary Smallest split index.
#' @param s0 Number of largest absolute coordinates in each norm.
#' @param powers Positive norm powers, including Inf.
#' @param seed Optional seed; NULL uses the caller's random stream.
#' @return List with p.value, minimum-p selected location and component results.
#' @details This comparator uses a multiplier bootstrap. HDLCP's proposed tests
#' do not bootstrap observations. This function uses the comparator's
#' finite-B p-value convention, including possible zero values.
#' @export
lzzl20_test <- function(X,B=500L,boundary=floor(.2*nrow(X)),s0=floor(ncol(X)/2),
                         powers=c(1:5,Inf),seed=NULL) {
  X <- .data(X); B <- .integer(B,"B",2); boundary <- .integer(boundary,"boundary")
  s0 <- .integer(s0,"s0")
  if(any(is.na(powers)|powers<=0)||!length(powers)||boundary>nrow(X)/2||s0>ncol(X))stop("Invalid norm settings.")
  work <- function() {
    z <- lzzl20_boot_stats_cpp(X,matrix(rnorm(nrow(X)*B),nrow(X),B),boundary,s0,powers)
    cp <- rowSums(sweep(z$boot_stats,1,z$stat,">"))/(B+1)
    bootp <- t(apply(z$boot_stats,1,function(x)(B-rank(x,ties.method="average"))/B))
    pv <- sum(apply(bootp,2,min)<=min(cp))/(B+1)
    list(p.value=pv,location=z$th[which.min(cp)],component.p=setNames(cp,powers),statistic=z$stat)
  }
  if(is.null(seed))work() else .with_seed(seed,work)
}

.all_methods <- function(X,ref,competitors,bootstrap_B) {
  fit <- hdlcp_test(X,ref)$tests[,c("method","p.value","location")]
  if(competitors) {
    methods <- list(DMS0=dms_test(X,0),DMS0.5=dms_test(X,.5),
                    ZWS22=zws22_test(X),LZZL20=lzzl20_test(X,B=bootstrap_B))
    fit <- rbind(fit,do.call(rbind,lapply(names(methods),function(nm)
      data.frame(method=nm,p.value=methods[[nm]]$p.value,location=methods[[nm]]$location))))
  }
  rownames(fit)<-NULL;fit
}
