# Formula checks against v37, using fixed arrays and numerical integration only.
# No Gaussian reference generation, bootstrap, or Monte Carlo replication.
audit_X <- function() outer(1:40,1:8,function(i,j)sin(i*j/7)+cos(i/3+j))

audit_reference <- function(k=1L,q=.25) {
  tails<-c(0,.025,.05,.5,.975,1)
  structure(list(n=40,p=8,k=k,q=q,B=39,lambda=3,seed=1,floor_level=-5,
    version="0.1.2",theta_grid=c(0,.3),tails=tails,noise=.02,
    lookups=lapply(q,HDLCP:::build_q_lookup),
    normal=matrix(cos(seq_len(39*39)),39,39),
    fixed=setNames(lapply(k,function(kk)list(
      anchor=cbind(c(10,8,6,3,1,0),c(25,20,15,8,2,0)),
      tangent=array(0,c(6,2,2)),lookup=HDLCP:::build_fixed_cluster_lookup(
        kk,8,116/39^2,seq(0,.95,.005)))),as.character(k))),class="hdlcp_reference")
}

test_that("Gaussian marginal levels and their densities satisfy (17)-(18)", {
  kap<-596/199^2
  for(u in c(1/400,5/200,.999999)) {
    b<-HDLCP:::.fixed_level(u,kap)
    expect_equal(HDLCP:::fixed_marginal_tail(b,kap),u,tolerance=1e-9)
    if(u<.99) {
      h<-1e-5
      derivative<-(HDLCP:::fixed_marginal_tail(b-h,kap)-
        HDLCP:::fixed_marginal_tail(b+h,kap))/(2*h)
      expect_equal(HDLCP:::fixed_marginal_density(b,kap),derivative,tolerance=1e-8)
    }
  }
  expect_equal(HDLCP:::fixed_marginal_tail(c(0,Inf),kap),c(1,0))
  expect_equal(HDLCP:::fixed_marginal_density(Inf,kap),0)
})

test_that("Gamma and its derivative agree with independent adaptive integration", {
  for(q in c(.1,.5,.9)) {
    lu<-HDLCP:::build_q_lookup(q)
    a<-qnorm(q/2,lower.tail=FALSE);mu<-q+2*a*dnorm(a)
    for(rho in c(.2,.7,.99975)) {
      sd<-sqrt(1-rho^2)
      deriv<-function(x) {
        m<-rho*x
        8*x*dnorm(x)*(m*(pnorm((a-m)/sd,lower.tail=FALSE)+pnorm((-a-m)/sd))+
          sd*(dnorm((a-m)/sd)-dnorm((-a-m)/sd)))
      }
      direct<-integrate(deriv,a,Inf,rel.tol=1e-10)$value
      expect_equal(HDLCP:::gamma_prime_at(lu,rho),direct,tolerance=1e-6)
      expect_equal(HDLCP:::gamma_prime_at(lu,-rho),-direct,tolerance=1e-6)
    }
    expect_equal(HDLCP:::gamma_prime_at(lu,1),4*mu)
    expect_equal(HDLCP:::gamma_at(lu,1),lu$tau2)
    expect_equal(HDLCP:::gamma_at(lu,0),0)
    # Gamma(r)=integral_0^r Gamma'(u)du, independent outer integration.
    integrated<-integrate(function(r)sapply(r,function(rr) {
      sd<-sqrt(1-rr^2)
      integrate(function(x)8*x*dnorm(x)*(rr*x*(pnorm((a-rr*x)/sd,lower.tail=FALSE)+
        pnorm((-a-rr*x)/sd))+sd*(dnorm((a-rr*x)/sd)-dnorm((-a-rr*x)/sd))),
        a,Inf,rel.tol=1e-8)$value
    }),0,.5,rel.tol=1e-7)$value
    expect_equal(HDLCP:::gamma_at(lu,.5),integrated,tolerance=1e-7)
  }
})

test_that("strong/weak rank clustering equals (19) with raw MOM mass", {
  p<-20;R<-diag(p);R[1,2]<-R[2,1]<-1.1;R[3,4]<-R[4,3]<-.6
  sm<-list(R=R,s2=.3,threshold=.3)
  lookup<-HDLCP:::build_fixed_cluster_lookup(3,p,596/199^2,seq(0,.95,.005))
  ref<-list(p=p,noise=.02,fixed=list(`3`=list(lookup=lookup)))
  direct<-mean(vapply(lookup,function(x) {
    u<-x$rank/p;b<-x$threshold;uG<-pchisq(b,1,lower.tail=FALSE)
    pair<-vapply(c(.95,.6),function(r) {
      a<-sqrt(b);sd<-sqrt(1-r*r)
      2*integrate(function(z)dnorm(z)*(pnorm((r*z-a)/sd)+
        pnorm((-r*z-a)/sd)),a,Inf,rel.tol=1e-9)$value
    },1)
    strong<-sum(pmax(pmin(pair*(u/uG)^2,u)-u*u,0))
    weakmass<-max(p*.3/2-sum(pmax(c(1.1,.6)^2-.02,0)),0)
    weak<-2*b*dnorm(sqrt(b))^2*(u/uG)^2*weakmass
    A<-2*(strong+weak)/(p*u)
    min(.3,A/(2+A))
  },1))
  expect_equal(HDLCP:::.theta(ref,3,sm),direct,tolerance=1e-8)
})

test_that("full kernel, variance, and crossing intensity use coherent corrections", {
  lu<-HDLCP:::build_q_lookup(.25);R<-diag(4)
  R[1,2]<-R[2,1]<--.6;R[3,4]<-R[4,3]<-.1
  sm<-list(R=R,s2=2*(.6^2+.1^2)/4,threshold=.2)
  t<-(1:24)/25;v<-t*(1-t);pr<-list(n=25,p=4,v=v)
  z<-HDLCP:::.kernel(pr,sm,lu)
  a<-(outer(t,t,pmin)-tcrossprod(t))/sqrt(outer(v,v))
  C<-HDLCP:::gamma_at(lu,a)+.5*(HDLCP:::gamma_at(lu,.6*a)+2*lu$mu^2*(.1*a)^2)
  omega<-lu$tau2+.5*(HDLCP:::gamma_at(lu,.6)+2*lu$mu^2*.1^2)
  D<-2*lu$mu+.25*(.6*HDLCP:::gamma_prime_at(lu,.6)+4*lu$mu^2*.1^2)
  expect_equal(z$omega,omega)
  expect_equal(z$eta,D/omega)
  expect_equal(tcrossprod(z$factor),C/omega,tolerance=1e-7)
  expect_equal(HDLCP:::.kernel(pr,sm,lu,FALSE)$omega,omega)
  sm$R<-abs(sm$R)
  expect_equal(HDLCP:::.kernel(pr,sm,lu,FALSE)$eta,z$eta)
})

test_that("Poisson paths retain returns and fixed geometric multiplicities", {
  paths<-cbind(c(3,-7,5,1),c(1,2,-6,4),c(-2,0,1,0))
  mult<-c(2L,1L,3L);orders<-c(1L,3L,5L)
  v<-c(.09,.21,.25,.16);trim<-2:3
  expanded<-paths[,rep(1:3,mult),drop=FALSE]
  direct<-unlist(lapply(orders,function(k) {
    sums<-apply(expanded,1,function(x)sum(sort(x,decreasing=TRUE)[seq_len(k)]))
    raw<-k*7+2*sums
    c(max(v*raw),max(raw[trim]))
  }))
  expect_equal(HDLCP:::tangent_path_stats_cpp(paths,mult,orders,7,2,v,trim),direct)
  # Zero draws checks the multi-order native interface without generating paths.
  empty<-HDLCP:::tangent_bank_cpp(0L,.2,orders,7,2,c(-2,-1,0,1),v,trim,-5)
  expect_equal(dim(empty),c(0L,6L))
  # Detailed balance for the two-sided Brownian marked Poisson measure (9).
  for(dt in c(.1,1,3))for(x in c(-2,.5))for(y in c(-1,2)) {
    forward<-exp(-x)*dnorm(y-x,mean=-dt/2,sd=sqrt(dt))
    backward<-exp(-y)*dnorm(x-y,mean=-dt/2,sd=sqrt(dt))
    expect_equal(forward,backward,tolerance=1e-12)
  }
})

test_that("critical inversion and Gaussian survival have the paper's endpoints", {
  crit<-c(5,4,3);tails<-c(0,.05,1)
  expect_equal(HDLCP:::.curve_p(6,crit,tails),0)
  expect_equal(HDLCP:::.curve_p(2,crit,tails),1)
  expect_equal(HDLCP:::.curve_p(3,crit,tails),1)
  expect_equal(HDLCP:::.curve_p(4,crit,tails),.05)
  expect_equal(HDLCP:::.curve_p(5,c(5,5,3),tails),.05)
  expect_equal(HDLCP:::.reference_tail(2,c(1,2,2,3)),.25)
  expect_equal(HDLCP:::.reference_tail(4,c(1,2,2,3)),0)
  expect_equal(HDLCP:::.reference_tail(0,c(1,2,2,3)),1)
  ref<-audit_reference(q=numeric());X<-audit_X()
  stat<-max(HDLCP:::.prepare(X,3)$sums[,1])
  ref$tails<-tails
  ref$fixed[[1]]$anchor<-cbind(c(stat+1,stat,stat-1),c(9,8,7))
  ref$fixed[[1]]$tangent<-array(0,c(3,2,2))
  z<-hdlcp_fixed_test(X,ref,alpha=.05)
  expect_equal(z$p.value,.05)
  expect_false(z$reject)
})

test_that("proportional scans, centering and skewness are (25)-(29)", {
  X<-audit_X();ref<-audit_reference(k=integer());q<-.25
  fit<-hdlcp_test(X,ref,keep_paths=TRUE)
  pr<-HDLCP:::.prepare(X,3);sm<-HDLCP:::.summaries(pr)
  lu<-ref$lookups[[1]];ker<-HDLCP:::.kernel(pr,sm,lu)
  kap<-116/39^2;Bq<-lu$mu+.5*lu$uq^2*dchisq(lu$uq,1)
  center<-8*lu$mu+(ceiling(8*q)-8*q)*lu$uq+8*kap*Bq
  U<-apply(pr$D,1,function(d)sum(pmax(d*d-lu$uq,0)+q*lu$uq-lu$mu)/sqrt(8))
  rho<-8^(-.5)+40^(-.5)+sqrt(log(320)/40)+log(320)^2/40
  skew<-sqrt(39*38)/37*mean((U-mean(U))^3)/(mean((U-mean(U))^2)+rho^2)^1.5
  skew<-min(rho,max(-rho,skew))
  expect_equal(unname(fit$diagnostics$q0.25['skewness']),skew)
  base<-(pr$sums[,2]-pr$v*center)/sqrt(8*ker$omega)
  expect_equal(fit$paths[['P0.25-0']],base)
  expect_equal(fit$paths[['P0.25-0.5']][3:37],(base/pr$v)[3:37])
  expect_true(all(is.na(fit$paths[['P0.25-0.5']][c(1,2,38,39)])))
  G<-ker$factor%*%ref$normal
  maxref<-apply(pr$v*(G+skew*(G^2-1)/6),2,max)
  expect_equal(fit$components$p.value[1],mean(maxref>max(base)))
  Z<-sqrt(2*log(pr$Tn))*max((base/pr$v)[3:37])-
    2*log(pr$Tn)-.5*log(log(pr$Tn))-log(ker$eta/sqrt(pi))
  expect_equal(fit$components$p.value[2],-expm1(-exp(-Z)))
  expect_equal(fit$components$location,c(which.max(base),2+which.max((base/pr$v)[3:37])))
})

test_that("Cauchy score and minimum-p location remain separate operations", {
  p<-c(.02,.17,.4,.7);w<-c(.1,.2,.3,.4)
  expect_equal(HDLCP:::.cauchy_result(p,w)$statistic,sum(w*tan(pi*(.5-p))))
  expect_equal(cauchy_combine(c(1e-310,.5),c(1e-308,1)),atan(pi*.01)/pi,tolerance=1e-10)
  expect_equal(HDLCP:::.gumbel_log_tail(c(40,1000)),c(-40,-1000))
  fit<-hdlcp_test(audit_X(),audit_reference(k=c(1,3),q=c(.25,.75)),keep_paths=TRUE)
  for(g in c(0,.5)) {
    co<-fit$components[fit$components$gamma==g,];cl<-fit$tests[fit$tests$method==paste0('CL',g),]
    expect_equal(cl$p.value,cauchy_combine(co$p.value))
    expect_equal(cl$statistic,HDLCP:::.cauchy_result(co$p.value)$statistic)
    expect_equal(cl$selected,co$method[which.min(co$log.p.value)])
    expect_equal(cl$location,co$location[which.min(co$log.p.value)])
  }
  old<-audit_reference();old$version<-'0.1.1'
  expect_error(hdlcp_test(audit_X(),old),'Rebuild')
})

test_that("table construction preserves batch, theta and order axes without random draws", {
  calls<-integer()
  testthat::local_mocked_bindings(
    .parallel=function(tasks,fun,workers,seed) {
      if(seed==11L)return(lapply(tasks,function(ids)outer(ids,1:4,'+')))
      lapply(tasks,fun)
    },
    .with_seed=function(seed,fun) {
      if(seed==13L)return(.02)
      if(seed==14L)return(matrix(cos(seq_len(24*39)),24,39))
      stop('Unexpected random reference request')
    },
    tangent_bank_cpp=function(B,theta,orders,beta,scale,s,v,trim,floor_level) {
      calls<<-c(calls,B)
      outer(seq_len(B),seq_len(2*length(orders)),'+')+theta
    },.package='HDLCP')
  ref<-hdlcp_calibrate(25,8,k=c(1,3),q=.25,B=39,workers=25,seed=11)
  expect_equal(sum(calls),39*8)
  expect_equal(sort(unique(calls)),c(7L,32L))
  expect_equal(range(ref$tails),c(0,1))
  expect_equal(dim(ref$normal),c(24L,39L))
  for(k in c(1,3))for(g in c(0,.5)) {
    base<-HDLCP:::.fixed_curve(ref,k,0,g)
    expect_equal(HDLCP:::.fixed_curve(ref,k,.075,g),base+.075)
    expect_equal(HDLCP:::.fixed_curve(ref,k,.125,g),base+.125)
  }
})
