#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

namespace {
void add_height(std::vector<double>& top,int time,int kmax,
                double height,int multiplicity) {
  const int start=time*kmax;
  for(int copy=0;copy<std::min(kmax,multiplicity);++copy) {
    if(height<=top[start+kmax-1])break;
    int h=kmax-1;
    while(h>0 && height>top[start+h-1]) {
      top[start+h]=top[start+h-1];--h;
    }
    top[start+h]=height;
  }
}

NumericVector scan_top(const std::vector<double>& top,IntegerVector orders,
                      double beta,double scale,NumericVector v,
                      const std::vector<bool>& keep_time) {
  const int kmax=*std::max_element(orders.begin(),orders.end());
  NumericVector out(2*orders.size(),R_NegInf);
  for(int i=0;i<v.size();++i) {
    std::vector<double> sums(kmax);double total=0;
    for(int h=0;h<kmax;++h){total+=top[i*kmax+h];sums[h]=total;}
    for(int a=0;a<orders.size();++a) {
      const double raw=orders[a]*beta+scale*sums[orders[a]-1];
      out[2*a]=std::max(out[2*a],v[i]*raw);
      if(keep_time[i])out[2*a+1]=std::max(out[2*a+1],raw);
    }
  }
  return out;
}
}

// Deterministic evaluation also used to check the path/rank formulas.
// [[Rcpp::export]]
NumericVector tangent_path_stats_cpp(NumericMatrix paths,IntegerVector mult,
 IntegerVector orders,double beta,double scale,NumericVector v,IntegerVector trim) {
 const int nt=paths.nrow(),kmax=*std::max_element(orders.begin(),orders.end());
 if(paths.ncol()!=mult.size() || nt!=v.size())stop("Invalid tangent paths.");
 std::vector<bool> keep(nt,false);
 for(int i:trim)keep.at(i-1)=true;
 std::vector<double> top(nt*kmax,R_NegInf);
 for(int j=0;j<paths.ncol();++j)
   for(int i=0;i<nt;++i)add_height(top,i,kmax,paths(i,j),mult[j]);
 return scan_top(top,orders,beta,scale,v,keep);
}

// [[Rcpp::export]]
NumericMatrix tangent_bank_cpp(int B,double theta,IntegerVector orders,double beta,
 double scale,NumericVector s,NumericVector v,IntegerVector trim,double floor_level) {
 const int nt=s.size(),kmax=*std::max_element(orders.begin(),orders.end());
 if(nt!=v.size() || nt<2 || theta<0 || theta>=1)stop("Invalid tangent grid.");
 std::vector<bool> keep(nt,false);
 for(int i:trim)keep.at(i-1)=true;
 std::vector<double> dt(nt-1),sd(nt-1);
 for(int i=0;i<nt-1;++i) {
   dt[i]=s[i+1]-s[i];
   if(dt[i]<=0)stop("Tangent times must increase.");
   sd[i]=std::sqrt(dt[i]);
 }
 const double intensity=(1-theta)*std::exp(-floor_level);
 NumericMatrix out(B,2*orders.size());
 for(int b=0;b<B;++b) {
   checkUserInterrupt();
   std::vector<double> top(nt*kmax,R_NegInf),path(nt);
   // Disjoint first-appearance roots partition the same Poisson path measure.
   // At root r, sample centers above the floor and retain only paths below it
   // at every earlier grid time. The stationary intensity exp(-x) is reversible
   // for a Brownian increment N(-dt/2,dt), in both time directions.
   for(int r=0;r<nt;++r) {
     const int candidates=R::rpois(intensity);
     for(int j=0;j<candidates;++j) {
       path[r]=floor_level+R::rexp(1);
       bool first=true;
       for(int i=r-1;i>=0;--i) {
         path[i]=path[i+1]+R::rnorm(-dt[i]/2,sd[i]);
         if(path[i]>=floor_level){first=false;break;}
       }
       if(!first)continue;
       const int multiplicity=1+R::rgeom(1-theta);
       // Retain below-floor excursions and later returns of the same center.
       for(int i=r+1;i<nt;++i)
         path[i]=path[i-1]+R::rnorm(-dt[i-1]/2,sd[i-1]);
       for(int i=0;i<nt;++i)add_height(top,i,kmax,path[i],multiplicity);
     }
   }
   bool sufficient=true;
   for(int i=0;i<nt;++i)
     if(top[i*kmax+kmax-1]<floor_level){sufficient=false;break;}
   if(!sufficient){for(int a=0;a<out.ncol();++a)out(b,a)=NA_REAL;continue;}
   NumericVector value=scan_top(top,orders,beta,scale,v,keep);
   for(int a=0;a<out.ncol();++a)out(b,a)=value[a];
 }
 return out;
}
