.root <- function(p,covariance,rho) {
  i <- seq_len(p)
  S <- switch(covariance,ar=rho^abs(outer(i,i,"-")),
    block=ifelse(outer(ceiling(i/5),ceiling(i/5),"=="),rho,0),independent=diag(p))
  diag(S)<-1
  e <- eigen(S,symmetric=TRUE)
  sweep(e$vectors,2,sqrt(pmax(e$values,0)),"*")%*%t(e$vectors)
}

.generate <- function(n,p,root,distribution,s,tau,delta) {
  Z <- if(distribution=="normal")matrix(rnorm(n*p),n,p) else matrix((rgamma(n*p,shape=4)-4)/2,n,p)
  X <- Z%*%root
  if(s>0) {
    support <- unique(as.integer(round(seq(1,p,length.out=s))))
    X[(floor(n*tau)+1):n,support] <- X[(floor(n*tau)+1):n,support,drop=FALSE]+delta
  }
  X
}

#' Generate a manuscript independent-component design
#' @param n,p Dimensions.
#' @param distribution Either normal or standardized Gamma(4,1).
#' @param covariance AR, five-coordinate block, or independent.
#' @param rho Correlation, in [0,1).
#' @param s Number of changed coordinates, zero under the null.
#' @param tau Change fraction strictly inside (0,1).
#' @param delta Per-coordinate shift after the change.
#' @param seed Optional integer seed.
#' @return Numeric n by p matrix. Uses the symmetric square root of Sigma,
#' which matters for non-Gaussian independent components.
#' @export
#' @examples
#' X <- simulate_hdlcp_data(40,20,s=3,delta=.8,seed=1)
simulate_hdlcp_data <- function(n=200L,p=200L,distribution=c("normal","gamma"),
  covariance=c("ar","block","independent"),rho=.5,s=0L,tau=.5,delta=0,seed=NULL) {
  n <- .integer(n,"n",25); p <- .integer(p,"p",2); s <- .integer(s,"s",0)
  distribution <- match.arg(distribution);covariance <- match.arg(covariance)
  if(length(rho)!=1||!is.finite(rho)||rho<0||rho>=1||s>p)stop("Invalid rho or s.")
  if(length(tau)!=1||!is.finite(tau)||tau<=0||tau>=1||floor(n*tau)<1||floor(n*tau)>=n)stop("Invalid tau.")
  if(length(delta)!=1||!is.finite(delta))stop("Invalid delta.")
  root <- .root(p,covariance,rho)
  fun <- function().generate(n,p,root,distribution,s,tau,delta)
  if(is.null(seed))fun() else .with_seed(seed,fun)
}

#' Paper simulation specifications
#' @return A list of all designs, replication counts, and signal grids in v37.
#' @export
paper_designs <- function() {
  list(n=200L,p=c(200L,400L),distribution=c("normal","gamma"),covariance=c("ar","block"),
    rho=.5,alpha=.05,k=c(1L,5L),q=c(.1,.9),gamma=c(0,.5),tau=c(.5,.3),
    off_center_multiplier=1.2,size_reps=1000L,order_null_reps=10000L,order_power_reps=10000L,
    comparison_null_reps=1000L,comparison_power_reps=1000L,location_reps=1000L,
    fixed=data.frame(s=1:10,delta=.65*(1:10)^(-.23)),
    proportional=data.frame(a=(1:9)/10,delta=c(.288,.207,.08*((3:9)/10)^(-.7))),
    comparison=data.frame(s=c(1,5,10,20,30,40,60,80,100,120,160,180),
      delta=c(.750,.430,.340,.270,.240,.210,.180,.150,.135,.120,.112,.112)),
    lzzl=list(B=500L,s0_fraction=.5,powers=c(1:5,Inf),boundary_fraction=.2))
}

.mc_methods <- function(n,p,root,distribution,s,tau,delta,reps,ref,competitors,B,workers,seed) {
  .parallel(as.list(seq_len(reps)),function(i) {
    X <- .generate(n,p,root,distribution,s,tau,delta)
    .all_methods(X,ref,competitors,B)
  },workers,seed)
}

.aggregate_methods <- function(values,cutoffs,truth=NULL,alpha=.05) {
  methods <- values[[1]]$method
  pv <- vapply(values,function(x)x$p.value,numeric(length(methods)))
  loc <- vapply(values,function(x)x$location,numeric(length(methods)))
  cutoff <- if(is.null(cutoffs))rep(alpha,length(methods)) else unname(cutoffs[methods])
  rejection <- rowMeans(pv<cutoff)
  err <- if(is.null(truth))matrix(NA_real_,nrow(pv),ncol(pv)) else abs(loc-truth)
  data.frame(method=methods,replications=ncol(pv),rejection=rejection,
    mcse=sqrt(rejection*(1-rejection)/ncol(pv)),mean_absolute_error=rowMeans(err),
    rmse=sqrt(rowMeans(err^2)),critical_p=cutoff,row.names=NULL)
}

#' Recompute the empirical size experiment (Table 1)
#' @param reps Null replications per design.
#' @param n Sample size.
#' @param p_values Dimensions to evaluate.
#' @param distributions Innovation laws.
#' @param covariances Covariance structures.
#' @param reference_B Reference-table draws.
#' @param competitors Include DMS0,DMS0.5,ZWS22,LZZL20.
#' @param bootstrap_B LZZL20 multiplier draws.
#' @param workers Parallel workers (default 25).
#' @param seed Reproducible seed.
#' @return Consolidated data frame, not replication-level dumps.
#' @export
reproduce_size <- function(reps=1000L,n=200L,p_values=c(200L,400L),
  distributions=c("normal","gamma"),covariances=c("ar","block"),
  reference_B=2000L,competitors=TRUE,bootstrap_B=500L,workers=25L,seed=10101L) {
  reps <- .integer(reps,"reps")
  if(competitors).dms_required()
  if(any(!distributions%in%c("normal","gamma"))||any(!covariances%in%c("ar","block","independent")))stop("Invalid design.")
  out<-list();id<-0L
  for(p in p_values) {
    ref<-hdlcp_calibrate(n,p,B=reference_B,workers=workers,seed=seed+p)
    for(cv in covariances)for(dist in distributions) {
      id<-id+1L;root<-.root(p,cv,.5)
      z<-.mc_methods(n,p,root,dist,0,.5,0,reps,ref,competitors,bootstrap_B,workers,seed+1000L+id)
      out[[id]]<-cbind(n=n,p=p,covariance=cv,distribution=dist,
                       .aggregate_methods(z,NULL))
    }
  }
  do.call(rbind,out)
}

.order_scores <- function(X,family,orders,lookups) {
  n<-nrow(X);p<-ncol(X);pr<-.prepare(X,ceiling(n^.3))
  if(family=="fixed")return(unlist(lapply(orders,function(k)
    c(max(pr$sums[,k]),max(pr$sums[pr$trim,k]/pr$v[pr$trim])))))
  sm<-.summaries(pr);kap<-(3*(n-1)-1)/(n-1)^2
  unlist(lapply(seq_along(orders),function(i) {
    q<-orders[i];lu<-lookups[[i]];kp<-ceiling(q*p)
    rho<-pmin(abs(sm$R[upper.tri(sm$R)]),1);rho<-rho[rho>sm$threshold]
    omega<-lu$tau2+2*lu$mu^2*sm$s2
    full<-omega+2/p*sum(gamma_at(lu,rho)-2*lu$mu^2*rho^2)
    slope<-2*lu$mu+2*lu$mu^2*sm$s2+1/p*sum(rho*gamma_prime_at(lu,rho)-4*lu$mu^2*rho^2)
    bias<-if(q==1)0 else .5*lu$uq^2*dchisq(lu$uq,1)
    center<-p*lu$mu+(kp-p*q)*lu$uq+p*kap*(lu$mu+bias)
    z<-(pr$sums[,kp]-pr$v*center)
    th<-max((z/(pr$v*sqrt(p*full)))[pr$trim])
    c(max(z/sqrt(p*full)),sqrt(2*log(pr$Tn))*th-
      (2*log(pr$Tn)+.5*log(log(pr$Tn))+log((slope/full)/sqrt(pi))))
  }))
}

#' Recompute order-selection power curves (Figures 1 and 2)
#' @param family Fixed k or proportional q.
#' @param reps Alternative replications per sparsity and location.
#' @param null_reps Independent H0 replications used once for critical values.
#' @param n,p Dimensions.
#' @param tau Change fractions.
#' @param workers Parallel workers.
#' @param seed Seed.
#' @return List with summary, critical values, and settings. Includes paired
#' Monte Carlo standard errors for power gains.
#' @details q=1 is calibrated directly by empirical H0 quantiles; no proportional
#' limit theorem is used at this endpoint. Results are fresh Monte Carlo draws.
#' @export
reproduce_order_power <- function(family=c("fixed","proportional"),reps=10000L,
  null_reps=10000L,n=200L,p=200L,tau=c(.5,.3),workers=25L,seed=20201L) {
  family<-match.arg(family)
  reps<-.integer(reps,"reps");null_reps<-.integer(null_reps,"null_reps")
  orders<-if(family=="fixed")1:10 else (1:10)/10
  if(p<10)stop("Order experiment needs p >= 10.")
  lookups<-if(family=="proportional")lapply(orders,build_q_lookup) else NULL
  root<-.root(p,"ar",.5)
  null<-.parallel(as.list(seq_len(null_reps)),function(i)
    .order_scores(.generate(n,p,root,"normal",0,.5,0),family,orders,lookups),workers,seed)
  critical<-apply(do.call(rbind,null),2,quantile,probs=.95,type=8)
  design<-if(family=="fixed")paper_designs()$fixed else paper_designs()$proportional
  if(family=="proportional")design$s<-round(p*design$a)
  out<-list();id<-0L;baseline<-if(family=="fixed")1L else length(orders)
  for(tt in tau)for(ii in seq_len(nrow(design))) {
    id<-id+1L;ss<-design$s[ii];delta<-design$delta[ii]*if(tt==.3)1.2 else 1
    z<-.parallel(as.list(seq_len(reps)),function(i) {
      score<-.order_scores(.generate(n,p,root,"normal",ss,tt,delta),family,orders,lookups)
      as.numeric(score>critical)
    },workers,seed+id)
    Z<-do.call(rbind,z);power<-colMeans(Z)
    basecols<-2*baseline-c(1,0)
    gains<-Z-Z[,rep(basecols,length(orders)),drop=FALSE]
    out[[id]]<-data.frame(family=family,tau=tt,s=ss,a=ss/p,order=rep(orders,each=2),
      gamma=rep(c(0,.5),length(orders)),delta=delta,power=power,gain=colMeans(gains),
      mcse=sqrt(power*(1-power)/reps),gain_mcse=apply(gains,2,sd)/sqrt(reps),replications=reps)
  }
  list(summary=do.call(rbind,out),critical=critical,
       settings=list(n=n,p=p,null_reps=null_reps,seed=seed))
}

#' Recompute adaptive power and localization (Figures 3 and 4)
#' @param reps Alternative replications per design.
#' @param null_reps Independent null replications for size-corrected power.
#' @param n,p Dimensions.
#' @param tau Change fractions.
#' @param sparsity Changed coordinate counts; paper uses twelve points.
#' @param delta Matching baseline signal vector.
#' @param reference_B Reference draws.
#' @param competitors Include all four competitors.
#' @param bootstrap_B LZZL20 multiplier draws.
#' @param workers Worker count.
#' @param seed Seed.
#' @return List with power/location summary, H0 cutoffs, and settings.
#' @export
reproduce_comparison <- function(reps=1000L,null_reps=1000L,n=200L,p=200L,tau=c(.5,.3),
  sparsity=paper_designs()$comparison$s,delta=paper_designs()$comparison$delta,
  reference_B=2000L,competitors=TRUE,bootstrap_B=500L,workers=25L,seed=30301L) {
  reps<-.integer(reps,"reps");null_reps<-.integer(null_reps,"null_reps")
  if(length(sparsity)!=length(delta)||any(sparsity>p)||any(sparsity<1))stop("Invalid signal grid.")
  if(competitors).dms_required()
  root<-.root(p,"ar",.5)
  ref<-hdlcp_calibrate(n,p,B=reference_B,workers=workers,seed=seed)
  null<-.mc_methods(n,p,root,"normal",0,.5,0,null_reps,ref,competitors,bootstrap_B,workers,seed+1L)
  cut<-setNames(apply(vapply(null,function(z)z$p.value,numeric(nrow(null[[1]]))),1,
    quantile,probs=.05,type=8),null[[1]]$method)
  out<-list();id<-0L
  for(tt in tau)for(i in seq_along(sparsity)) {
    id<-id+1L;dd<-delta[i]*if(tt==.3)1.2 else 1
    z<-.mc_methods(n,p,root,"normal",sparsity[i],tt,dd,reps,ref,competitors,bootstrap_B,workers,seed+100+id)
    out[[id]]<-cbind(tau=tt,s=sparsity[i],delta=dd,
      .aggregate_methods(z,cut,truth=floor(n*tt)))
  }
  list(summary=do.call(rbind,out),critical_p=cut,
       settings=list(n=n,p=p,null_reps=null_reps,seed=seed))
}

#' Plot newly computed paper simulation results
#' @param result Object returned by reproduce_order_power or reproduce_comparison.
#' @param file Output PDF path.
#' @param metric Power or localization for the comparison experiment.
#' @return Invisibly, file. Creates only the requested final figure.
#' @export
plot_paper_results <- function(result,file,metric=c("power","location")) {
  metric<-match.arg(metric);z<-result$summary
  pdf(file,width=12,height=if("gain"%in%names(z))7 else 5.5)
  on.exit(dev.off())
  for(tt in unique(z$tau)) {
    d<-z[z$tau==tt,]
    if("gain"%in%names(z)) {
      par(mfrow=if(d$family[1]=="fixed")c(2,5) else c(3,3),mar=c(3,3,2,1),oma=c(1,2,2,0))
      for(s in unique(d$s)) {
        a<-d[d$s==s,];o<-sort(unique(a$order))
        matplot(o,cbind(a$gain[a$gamma==0],a$gain[a$gamma==.5]),type="b",pch=c(16,17),
          lty=c(1,2),col=c("#0072B2","#D55E00"),xlab="Selected order",ylab="",main=paste("s =",s))
        abline(h=0,lty=3)
        if(s==unique(d$s)[1])legend("topright",c("gamma=0","gamma=0.5"),col=c("#0072B2","#D55E00"),lty=1:2,bty="n",cex=.7)
      }
      mtext(if(d$family[1]=="fixed")"Power(k)-Power(1)" else "Power(q)-Power(1)",2,outer=TRUE,line=.5)
      mtext(paste("tau =",tt),3,outer=TRUE)
    } else {
      par(mfrow=c(1,1),mar=c(4,4,3,1))
      methods<-intersect(c("CL0","CL0.5","DMS0","DMS0.5","ZWS22","LZZL20"),unique(d$method))
      ss<-unique(d$s);key<-if(metric=="power")"rejection" else "mean_absolute_error"
      vals<-sapply(methods,function(m)d[d$method==m,key])
      matplot(seq_along(ss),vals,type="b",xaxt="n",pch=seq_along(methods),lty=seq_along(methods),
        xlab="Number of changed coordinates s",ylab=if(metric=="power")"Size-corrected power" else "Mean absolute location error",
        main=paste("tau =",tt))
      axis(1,at=seq_along(ss),labels=ss);legend("topleft",methods,col=seq_along(methods),lty=seq_along(methods),bty="n",cex=.8)
    }
  }
  invisible(file)
}

#' Recompute all paper simulation designs
#' @param output_dir Directory for final consolidated summaries and figures only.
#' @param mode Paper replication counts or a small smoke run.
#' @param workers Worker count.
#' @param seed Seed.
#' @param competitors Include competitors.
#' @param reference_B Reference draws (paper default 2000; increase for precision).
#' @return Named list of final summaries. Never writes per-replication data or checkpoints.
#' @export
reproduce_paper <- function(output_dir="hdlcp-results",mode=c("paper","smoke"),workers=25L,
                            seed=40401L,competitors=TRUE,reference_B=NULL) {
  mode<-match.arg(mode);small<-mode=="smoke"
  reps<-if(small)2L else 1000L;ord<-if(small)2L else 10000L
  B<-if(is.null(reference_B))if(small)39L else 2000L else reference_B
  result<-list(size=reproduce_size(reps=reps,reference_B=B,competitors=competitors,
    bootstrap_B=if(small)9L else 500L,workers=workers,seed=seed),
    fixed=reproduce_order_power("fixed",reps=ord,null_reps=ord,workers=workers,seed=seed+1L),
    proportional=reproduce_order_power("proportional",reps=ord,null_reps=ord,workers=workers,seed=seed+2L),
    comparison=reproduce_comparison(reps=reps,null_reps=reps,reference_B=B,competitors=competitors,
      bootstrap_B=if(small)9L else 500L,workers=workers,seed=seed+3L))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  write.csv(result$size,file.path(output_dir,"table1_size.csv"),row.names=FALSE)
  for(nm in c("fixed","proportional","comparison"))
    write.csv(result[[nm]]$summary,file.path(output_dir,paste0(nm,"_summary.csv")),row.names=FALSE)
  plot_paper_results(result$fixed,file.path(output_dir,"figure1.pdf"))
  plot_paper_results(result$proportional,file.path(output_dir,"figure2.pdf"))
  plot_paper_results(result$comparison,file.path(output_dir,"figure3.pdf"))
  plot_paper_results(result$comparison,file.path(output_dir,"figure4.pdf"),"location")
  result
}
