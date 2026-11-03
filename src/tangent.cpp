#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix tangent_bank_cpp(int B, double theta, int k, double beta,
 double scale, NumericVector s, NumericVector v, IntegerVector trim, double floor_level) {
 NumericMatrix out(B,2); std::vector<bool> keep_time(s.size(),false);
 for(int i:trim) keep_time[i-1]=true;
 double intensity=(1-theta)*std::exp(-floor_level);
 for(int b=0;b<B;b++) {
   std::vector<double> x; std::vector<int> mult;
   int nc=R::rpois(intensity);
   for(int j=0;j<nc;j++){x.push_back(floor_level+R::rexp(1));mult.push_back(1+R::rgeom(1-theta));}
   double best0=R_NegInf,besth=R_NegInf;
   for(int i=0;i<s.size();i++) {
     std::vector<double> vals;
     for(size_t j=0;j<x.size();j++) for(int h=0;h<std::min(k,mult[j]);h++) vals.push_back(x[j]);
     if(vals.size()>=size_t(k)) {
       std::partial_sort(vals.begin(),vals.begin()+k,vals.end(),std::greater<double>());
       double sk=0;for(int h=0;h<k;h++)sk+=vals[h];
       double raw=k*beta+scale*sk; best0=std::max(best0,v[i]*raw);
       if(keep_time[i])besth=std::max(besth,raw);
     }
     if(i+1<s.size()) {
       double dt=s[i+1]-s[i];std::vector<double> nx;std::vector<int> nm;
       for(size_t j=0;j<x.size();j++){double z=x[j]+R::rnorm(-dt/2,std::sqrt(dt));
         if(z>=floor_level){nx.push_back(z);nm.push_back(mult[j]);}}
       int np=R::rpois(intensity);
       for(int j=0;j<np;j++){double z=floor_level+R::rexp(1);
         if(R::runif(0,1)<R::pnorm((floor_level-z+dt/2)/std::sqrt(dt),0,1,1,0)){
           nx.push_back(z);nm.push_back(1+R::rgeom(1-theta));}}
       x.swap(nx);mult.swap(nm);
     }
   }
   out(b,0)=best0;out(b,1)=besth;
 }
 return out;
}
