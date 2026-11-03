#' Read and preprocess the SECOM window used in the paper
#' @param data_file Path to UCI secom.data.
#' @param labels_file Path to UCI secom_labels.data.
#' @param rows Contiguous window in timestamp order, default 493:892.
#' @param missing_threshold Maximum global missing fraction.
#' @param flat_run Remove variables with runs at least this long within the window.
#' @return List with X, timestamps, labels and original feature indices.
#' @details Missingness screening, imputation and constant screening use the whole
#' timestamp-ordered data. Window-specific difference variance and flat-run
#' screening follow. No labels are used to select variables.
#' @export
read_secom <- function(data_file,labels_file,rows=493:892,missing_threshold=.01,flat_run=20L) {
  X<-as.matrix(read.table(data_file,na.strings=c("NaN","NA"),fill=TRUE))
  lines<-readLines(labels_file,warn=FALSE)
  labels<-as.integer(sub("[[:space:]].*$","",trimws(lines)))
  stamp<-gsub('"','',sub("^[^[:space:]]+[[:space:]]+","",trimws(lines)),fixed=TRUE)
  times<-as.POSIXct(stamp,format="%d/%m/%Y %H:%M:%S",tz="UTC")
  if(anyNA(times)||length(times)!=nrow(X))stop("Invalid SECOM timestamps or row count.")
  idx<-order(times);X<-X[idx,,drop=FALSE];times<-times[idx];labels<-labels[idx]
  keep<-which(colMeans(!is.finite(X))<=missing_threshold);X<-X[,keep,drop=FALSE]
  med<-apply(X,2,median,na.rm=TRUE);missing<-which(!is.finite(X),arr.ind=TRUE)
  if(nrow(missing))X[missing]<-med[missing[,2]]
  ok<-apply(X,2,var)>.Machine$double.eps & colMeans(diff(X)^2)>.Machine$double.eps
  X<-X[,ok,drop=FALSE];keep<-keep[ok]
  if(!length(rows)||any(rows<1|rows>nrow(X))||any(diff(rows)!=1))stop("Invalid contiguous rows.")
  X<-X[rows,,drop=FALSE]
  ok<-colMeans(diff(X)^2)>.Machine$double.eps
  X<-X[,ok,drop=FALSE];keep<-keep[ok]
  ok<-vapply(seq_len(ncol(X)),function(j)max(rle(X[,j])$lengths)<flat_run,logical(1))
  list(X=X[,ok,drop=FALSE],timestamps=times[rows],labels=labels[rows],
       features=keep[ok],flat_removed=keep[!ok],rows=rows)
}

#' Recompute the SECOM application (Tables 4-5 and Figures 5-6)
#' @param data_file,labels_file UCI data and label paths.
#' @param output_dir Directory for final summary tables and figures.
#' @param reference_B Gaussian reference draws.
#' @param reps Side-wise subsampling replications; zero skips Table 5.
#' @param side_sizes Samples drawn without replacement per side.
#' @param competitors Include DMS0,DMS0.5,ZWS22,LZZL20.
#' @param bootstrap_B LZZL20 multiplier draws.
#' @param workers Parallel workers.
#' @param seed Seed.
#' @return Full-window summary, Welch diagnostics, subsampling summary, preprocessing metadata.
#' @export
reproduce_secom <- function(data_file,labels_file,output_dir="hdlcp-secom",
  reference_B=2000L,reps=500L,side_sizes=c(50L,100L,150L),competitors=TRUE,
  bootstrap_B=500L,workers=25L,seed=50501L) {
  reps<-.integer(reps,"reps",0)
  if(competitors).dms_required()
  data<-read_secom(data_file,labels_file);X<-data$X;p<-ncol(X);truth<-229L
  if(nrow(X)!=400||p!=361)stop("SECOM preprocessing did not produce the paper's 400 by 361 window.")
  ref<-hdlcp_calibrate(400,p,B=reference_B,workers=workers,seed=seed)
  full<-.with_seed(seed+1,function().all_methods(X,ref,competitors,bootstrap_B))
  full$absolute_error<-abs(full$location-truth)
  welch<-do.call(rbind,lapply(seq_len(p),function(j) {
    z<-t.test(X[1:truth,j],X[(truth+1):400,j])
    data.frame(feature=data$features[j],t=unname(z$statistic),p.value=z$p.value)
  }))
  subs<-list()
  if(reps>0)for(m in side_sizes) {
    m<-.integer(m,"side size",13)
    if(m>171)stop("side size must not exceed 171.")
    rf<-hdlcp_calibrate(2*m,p,B=reference_B,workers=workers,seed=seed+m)
    z<-.parallel(as.list(seq_len(reps)),function(i) {
      rows<-c(sort(sample.int(truth,m)),truth+sort(sample.int(400-truth,m)))
      .all_methods(X[rows,,drop=FALSE],rf,competitors,bootstrap_B)
    },workers,seed+1000+m)
    subs[[as.character(m)]]<-cbind(m=m,n=2*m,.aggregate_methods(z,NULL,truth=m))
  }
  summary<-if(length(subs))do.call(rbind,subs) else data.frame()
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  write.csv(full,file.path(output_dir,"table4_secom.csv"),row.names=FALSE)
  write.csv(welch,file.path(output_dir,"welch_summary.csv"),row.names=FALSE)
  if(nrow(summary))write.csv(summary,file.path(output_dir,"table5_secom.csv"),row.names=FALSE)
  local({pdf(file.path(output_dir,"figure5.pdf"),width=10,height=8);on.exit(dev.off())
    par(mfrow=c(3,1),mar=c(3,4,2,1))
    for(j in head(order(abs(welch$t),decreasing=TRUE),3)) {
      plot(seq_len(400),X[,j],type="l",xlab="Ordered observation",ylab="Measurement",
           main=paste("Original feature",data$features[j]));abline(v=truth,lty=2)
      segments(1,mean(X[1:truth,j]),truth,mean(X[1:truth,j]),col=2,lwd=2)
      segments(truth+1,mean(X[(truth+1):400,j]),400,mean(X[(truth+1):400,j]),col=4,lwd=2)
    }})
  local({pdf(file.path(output_dir,"figure6.pdf"),width=10,height=4);on.exit(dev.off())
    par(mfrow=c(1,2));hist(welch$t,main="Welch t-statistics",xlab="t")
    hist(welch$p.value,main="Two-sided p-values",xlab="p-value")})
  list(full=full,welch=welch,subsampling=summary,features=data$features,
       flat_removed=data$flat_removed,timestamps=range(data$timestamps))
}
