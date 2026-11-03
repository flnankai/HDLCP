#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

namespace {

double safe_sqrt(double x) {
  return std::sqrt(std::max(x, 1e-12));
}

std::vector<double> top_values(const arma::rowvec& x, int s0) {
  const int p = x.n_elem;
  std::vector<double> vals(p);
  for (int j = 0; j < p; ++j) vals[j] = std::abs(x(j));
  if (s0 < p) {
    std::nth_element(vals.begin(), vals.begin() + s0, vals.end(),
                     std::greater<double>());
  }
  return vals;
}

void update_norms(const arma::rowvec& values, const NumericVector& powers,
                  int s0, int location, NumericVector& stat,
                  IntegerVector& th) {
  std::vector<double> vals = top_values(values, s0);
  const int n_power = powers.size();
  for (int ip = 0; ip < n_power; ++ip) {
    double norm_value = 0.0;
    if (R_finite(powers[ip])) {
      const double power = powers[ip];
      double sum_power = 0.0;
      for (int j = 0; j < s0; ++j) sum_power += std::pow(vals[j], power);
      norm_value = std::pow(sum_power, 1.0 / power);
    } else {
      for (int j = 0; j < s0; ++j) norm_value = std::max(norm_value, vals[j]);
    }
    if (norm_value > stat[ip]) {
      stat[ip] = norm_value;
      th[ip] = location;
    }
  }
}

void update_boot_norms(const arma::rowvec& values, const NumericVector& powers,
                       int s0, arma::vec& stat) {
  std::vector<double> vals = top_values(values, s0);
  const int n_power = powers.size();
  for (int ip = 0; ip < n_power; ++ip) {
    double norm_value = 0.0;
    if (R_finite(powers[ip])) {
      const double power = powers[ip];
      double sum_power = 0.0;
      for (int j = 0; j < s0; ++j) sum_power += std::pow(vals[j], power);
      norm_value = std::pow(sum_power, 1.0 / power);
    } else {
      for (int j = 0; j < s0; ++j) norm_value = std::max(norm_value, vals[j]);
    }
    if (norm_value > stat(ip)) stat(ip) = norm_value;
  }
}

}  // namespace

// [[Rcpp::export]]
List lzzl20_boot_stats_cpp(const arma::mat& X, const arma::mat& e,
                           int n_boundary, int s0,
                           NumericVector powers) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const int B = e.n_cols;
  const int n_power = powers.size();
  if (n_boundary < 1 || n_boundary > n / 2) {
    stop("n_boundary must be in [1, floor(n / 2)].");
  }
  if (s0 < 1 || s0 > p) stop("s0 must be in [1, p].");
  if (static_cast<int>(e.n_rows) != n) stop("The multiplier matrix has incompatible n.");

  arma::rowvec sum_x = arma::sum(X, 0);
  arma::rowvec sum_x2 = arma::sum(arma::square(X), 0);
  arma::rowvec s2_0 = sum_x2 - arma::square(sum_x) / double(n);

  arma::mat cum_x(n - 1, p, arma::fill::zeros);
  arma::mat sig(n - 1, p, arma::fill::zeros);
  arma::rowvec running_x(p, arma::fill::zeros);
  NumericVector observed(n_power);
  IntegerVector th(n_power);

  for (int t = 1; t < n; ++t) {
    running_x += X.row(t - 1);
    cum_x.row(t - 1) = running_x;
    const double tt = double(t) / double(n);
    const double v = tt * (1.0 - tt);
    arma::rowvec cusum0 = (running_x - double(t) * sum_x / double(n)) /
      std::sqrt(double(n));
    arma::rowvec cusum05 = cusum0 / std::sqrt(v);
    arma::rowvec sig_row(p, arma::fill::zeros);
    arma::rowvec standardized(p, arma::fill::zeros);
    for (int j = 0; j < p; ++j) {
      sig_row(j) = safe_sqrt((s2_0(j) - cusum05(j) * cusum05(j)) / double(n));
      standardized(j) = cusum0(j) / sig_row(j);
    }
    sig.row(t - 1) = sig_row;
    if (t >= n_boundary && t <= n - n_boundary) {
      update_norms(standardized, powers, s0, t, observed, th);
    }
  }

  arma::mat boot_stats(n_power, B, arma::fill::zeros);
  for (int b = 0; b < B; ++b) {
    arma::vec eb = e.col(b);
    const double total_e = arma::sum(eb);
    arma::rowvec total_ex = eb.t() * X;
    double left_e = 0.0;
    arma::rowvec left_ex(p, arma::fill::zeros);
    arma::vec current(n_power, arma::fill::zeros);

    for (int t = 1; t < n; ++t) {
      const double et = eb(t - 1);
      left_e += et;
      left_ex += et * X.row(t - 1);
      const int nl = t;
      const int nr = n - t;
      arma::rowvec left_x = cum_x.row(t - 1);
      arma::rowvec right_x = sum_x - left_x;
      arma::rowvec right_ex = total_ex - left_ex;
      const double right_e = total_e - left_e;
      arma::rowvec dl = left_ex / double(nl) -
        left_e * left_x / double(nl * nl);
      arma::rowvec dr = right_ex / double(nr) -
        right_e * right_x / double(nr * nr);
      const double coef = double(nl * nr) / (double(n) * std::sqrt(double(n)));
      arma::rowvec standardized = coef * (dl - dr);
      for (int j = 0; j < p; ++j) standardized(j) /= sig(t - 1, j);
      if (t >= n_boundary && t <= n - n_boundary) {
        update_boot_norms(standardized, powers, s0, current);
      }
    }
    boot_stats.col(b) = current;
  }

  return List::create(
    _["stat"] = observed,
    _["th"] = th,
    _["boot_stats"] = boot_stats
  );
}
