# Deterministic identities only. No reference simulation or Monte Carlo experiment.
fixture <- function() outer(seq_len(40),seq_len(8),function(i,j) sin(i*j/7)+cos(i/3+j))

test_that("CUSUM and all adjacent differences equal their direct formulas", {
  X<-fixture();z<-HDLCP:::.prepare(X,3)
  expected<-colMeans((X[-1,]-X[-40,])^2)/2
  expect_equal(z$sigma2,expected)
  for(r in c(1,13,39)) {
    w<-(colSums(X[seq_len(r),,drop=FALSE])-r/40*colSums(X))/sqrt(40*expected)
    expect_equal(z$W[r,],w)
    expect_equal(z$sums[r,5],sum(sort(w^2,decreasing=TRUE)[1:5]))
  }
})

test_that("MOM blocks have the prescribed separation and trace exclusions", {
  for(count in c(5,20)) {
    blocks<-HDLCP:::.blocks(39,count)
    expect_length(blocks,count)
    expect_equal(vapply(blocks[-1],min,1)-vapply(head(blocks,-1),max,1),rep(2,count-1))
  }
  z<-HDLCP:::.prepare(fixture(),3);sm<-HDLCP:::.summaries(z)
  traces<-vapply(HDLCP:::.blocks(39,5),function(i) {
    pairs<-expand.grid(i=i,j=i);pairs<-pairs[abs(pairs$i-pairs$j)>1,]
    mean(vapply(seq_len(nrow(pairs)),function(h)sum(z$D[pairs$i[h],]*z$D[pairs$j[h],])^2,1))
  },1)
  expect_equal(sm$s2,max(median(traces)/8-1,0))
  perm<-c(8,2,5,1,7,3,4,6)
  rotated<-HDLCP:::.summaries(HDLCP:::.prepare(fixture()[,perm],3))
  expect_equal(rotated$s2,sm$s2)
  expect_equal(rotated$R,sm$R[perm,perm])
})

test_that("Cauchy formula and appendix endpoints hold", {
  p<-c(.02,.13,.27)
  expect_equal(cauchy_combine(p),.5-atan(mean(tan(pi*(.5-p))))/pi)
  expect_equal(cauchy_combine(.3),.3)
  expect_equal(cauchy_combine(c(.25,.75)),.5)
  expect_equal(cauchy_combine(c(.25,.75),c(1e308,1e308)),.5)
  expect_equal(cauchy_combine(.5),.5)
  expect_equal(cauchy_combine(c(0,1)),0)
  expect_equal(cauchy_combine(c(.1,1)),1)
  expect_error(cauchy_combine(c(NA,.2)))
})

test_that("normalization and exact q=1 moments agree", {
  lu<-HDLCP:::build_q_lookup(1)
  expect_equal(lu$mu,1);expect_equal(lu$tau2,2)
  expect_equal(HDLCP:::gamma_at(lu,c(0,.5,1)),c(0,.5,2))
  kap<-(3*39-1)/39^2;x<-c(1,3,10)
  expect_equal(HDLCP:::fixed_marginal_tail(x,kap),pchisq(x,1,lower.tail=FALSE)+kap*x*(x+1)*dchisq(x,1)/4)
})

test_that("max-norm correlation projection solves a deterministic SDP", {
  skip_if_not_installed("CVXR",minimum_version="1.9.0")
  K<-matrix(c(1,.9,.9,.9,1,-.9,.9,-.9,1),3,3)
  P<-HDLCP:::.project_kernel(K)
  expect_equal(diag(P),rep(1,3),tolerance=1e-6)
  expect_true(min(eigen(P,symmetric=TRUE)$values)>-1e-6)
  expect_equal(max(abs(P-K)),.4,tolerance=1e-5)
})

test_that("critical envelope is nonincreasing and inversion respects thresholds", {
  ref<-list(theta_grid=c(0,.3),fixed=list(`1`=list(anchor=cbind(c(5,4,3),c(5,4,3)),
    tangent=array(c(rep(0,6),c(0,2,0,0,2,0)),c(3,2,2)))))
  curve<-HDLCP:::.fixed_curve(ref,1,.3,0)
  expect_equal(curve,c(6,6,3))
  expect_true(all(diff(curve)<=0))
  expect_equal(HDLCP:::.curve_p(4.5,c(5,4,3),c(.01,.05,.1)),.03)
})

test_that("symmetric covariance root generates the target matrix", {
  r<-HDLCP:::.root(8,"ar",.5)
  expect_equal(r,t(r))
  expect_equal(r%*%r,.5^abs(outer(1:8,1:8,"-")),tolerance=1e-12)
  b<-HDLCP:::.root(8,"block",.5);S<-b%*%b
  expect_equal(diag(S),rep(1,8));expect_equal(S[1,5],.5);expect_equal(S[1,6],0,tolerance=1e-12)
})

test_that("degenerate scales use the paper's convention", {
  ref<-structure(list(n=40,p=8,k=1,q=.5,B=39,seed=1,lambda=3),class="hdlcp_reference")
  fit<-hdlcp_test(matrix(1,40,8),ref)
  expect_true(all(fit$tests$p.value==1))
  expect_equal(unique(fit$tests$location),c(1L,3L))
})

test_that("reproduction defaults and native interfaces are present without running simulations", {
  d<-paper_designs()
  expect_equal(d$fixed$delta,.65*(1:10)^(-.23))
  expect_equal(d$comparison$s,c(1,5,10,20,30,40,60,80,100,120,160,180))
  expect_equal(d$off_center_multiplier,1.2)
  expect_true(is.function(HDLCP:::tangent_bank_cpp))
  expect_true(is.function(HDLCP:::lzzl20_boot_stats_cpp))
  expect_true(is.function(HDLCP:::zws22_fast_stats))
})

test_that("the complete interface and localization agree on a deterministic fixture", {
  X<-fixture();lu<-HDLCP:::build_q_lookup(.25)
  tails<-c(.025,.05,.5,.975)
  ref<-structure(list(n=40,p=8,k=1,q=.25,B=39,lambda=3,seed=1,floor_level=-5,
    theta_grid=c(0,.3),tails=tails,noise=.02,lookups=list(lu),
    normal=matrix(cos(seq_len(39*39)),39,39),
    fixed=list(`1`=list(anchor=cbind(c(8,6,3,1),c(20,15,8,2)),
      tangent=array(0,c(4,2,2)),lookup=HDLCP:::build_fixed_cluster_lookup(1,8,116/39^2,seq(0,.95,.005))))),
    class="hdlcp_reference")
  fit<-hdlcp_test(X,ref,keep_paths=TRUE)
  expect_equal(nrow(fit$tests),6)
  expect_true(all(fit$tests$p.value>=0 & fit$tests$p.value<=1))
  for(i in seq_len(nrow(fit$components))) {
    z<-fit$components[i,];path<-fit$paths[[z$method]]
    allowed<-if(z$gamma==0)1:39 else 3:37
    expect_equal(z$location,allowed[which.max(path[allowed])])
  }
  for(g in c(0,.5)) {
    z<-fit$components[fit$components$gamma==g,]
    best<-z$location[which.min(z$p.value)]
    expect_equal(fit$tests$location[fit$tests$method==paste0("CL",g)],best)
  }
})
