#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

namespace {

inline double sign_pow_minus_one(const int exponent) {
  return (exponent % 2 == 0) ? 1.0 : -1.0;
}

List zws_one_fast(const arma::mat& X, const int q, const int truncation) {
  const int n = static_cast<int>(X.n_rows);
  const int p = static_cast<int>(X.n_cols);
  const int first_k = 2 * truncation;
  const int last_k = n - 2 * truncation;

  if (q < 2 || truncation < q || first_k > last_k) {
    stop("Invalid q/truncation or sample size for ZWS22");
  }

  arma::vec factorial(q, arma::fill::zeros);
  arma::mat choose(q, n, arma::fill::zeros);
  for (int order = 1; order <= q; ++order) {
    factorial(order - 1) = R::gammafn(order + 1.0);
    for (int sample_size = order; sample_size <= n; ++sample_size) {
      choose(order - 1, sample_size - 1) = R::choose(sample_size, order);
    }
  }

  const int n_k = last_k - first_k + 1;
  const int denom_width = n - 4 * truncation + 2;
  arma::rowvec num_sum(n_k, arma::fill::zeros);
  arma::mat denom(n_k, denom_width, arma::fill::zeros);
  arma::mat presum(q, n, arma::fill::zeros);
  arma::mat postsum(q, n, arma::fill::zeros);
  arma::vec temp_sum(q, arma::fill::zeros);

  for (int variable = 0; variable < p; ++variable) {
    const arma::vec x = X.col(variable);
    const arma::rowvec x_row = x.t();
    presum.zeros();
    postsum.zeros();
    presum.row(0) = arma::cumsum(x_row);
    postsum.row(0) = arma::cumsum(arma::reverse(x_row));

    // Elementary symmetric prefix and suffix sums.  The upper limit is kept
    // identical to the DMS package implementation for bit-level comparability
    // of all entries that enter the statistic.
    for (int order = 2; order <= q; ++order) {
      for (int i = order; i <= n + order - q; ++i) {
        presum(order - 1, i - 1) =
          presum(order - 1, i - 2) +
          presum(order - 2, i - 2) * x(i - 1);
        postsum(order - 1, i - 1) =
          postsum(order - 1, i - 2) +
          postsum(order - 2, i - 2) * x(n - i);
      }
    }

    for (int ik = 0; ik < n_k; ++ik) {
      const int k = first_k + ik;
      double numerator =
        factorial(q - 1) * presum(q - 1, k - 1) *
          choose(q - 1, n - k - 1) +
        sign_pow_minus_one(q) * factorial(q - 1) *
          postsum(q - 1, n - k - 1) * choose(q - 1, k - 1);

      for (int c = 1; c <= q - 1; ++c) {
        numerator +=
          sign_pow_minus_one(q - c) * factorial(c - 1) *
          factorial(q - c - 1) * presum(c - 1, k - 1) *
          choose(c - 1, n - k - q + c - 1) *
          choose(q - c - 1, k - c - 1) *
          postsum(q - c - 1, n - k - 1);
      }
      num_sum(ik) += numerator;

      // Left self-normalizer.  The original code repeatedly summed an
      // x-subvector here; the prefix-sum difference is algebraically exact
      // and removes the accidental O(p n^3) work.
      for (int t = truncation; t <= k - truncation; ++t) {
        double dq = presum(q - 1, t - 1) *
          choose(q - 1, k - t - 1) * factorial(q - 1);
        temp_sum(0) = presum(0, k - 1) - presum(0, t - 1);

        for (int c = 2; c <= q; ++c) {
          dq += temp_sum(c - 2) * presum(q - c, t - 1) *
            sign_pow_minus_one(c - 1) *
            choose(q - c, k - t - c) *
            choose(c - 2, t - q + c - 2) *
            factorial(c - 2) * factorial(q - c);
          temp_sum(c - 1) =
            presum(c - 1, k - 1) - presum(c - 1, t - 1);
          for (int d = 1; d <= c - 1; ++d) {
            temp_sum(c - 1) -= presum(d - 1, t - 1) *
              temp_sum(c - d - 1);
          }
        }
        denom(ik, t - truncation) +=
          dq + temp_sum(q - 1) * sign_pow_minus_one(q) *
          choose(q - 1, t - 1) * factorial(q - 1);
      }

      // Right self-normalizer, again using the already available prefix sum.
      for (int t = truncation; t <= n - truncation - k; ++t) {
        double dq = postsum(q - 1, t - 1) * sign_pow_minus_one(q) *
          choose(q - 1, n - t - k - 1) * factorial(q - 1);
        temp_sum(0) = presum(0, n - t - 1) - presum(0, k - 1);

        for (int c = 2; c <= q; ++c) {
          dq += temp_sum(c - 2) * postsum(q - c, t - 1) *
            sign_pow_minus_one(q - c + 1) *
            choose(q - c, n - t - k - c) *
            choose(c - 2, t - q + c - 2) *
            factorial(c - 2) * factorial(q - c);
          temp_sum(c - 1) =
            postsum(c - 1, n - k - 1) - postsum(c - 1, t - 1);
          for (int d = 1; d <= c - 1; ++d) {
            temp_sum(c - 1) -= postsum(d - 1, t - 1) *
              temp_sum(c - d - 1);
          }
        }
        denom(ik, t + k - 3 * truncation + 1) +=
          dq + temp_sum(q - 1) * choose(q - 1, t - 1) *
          factorial(q - 1);
      }
    }
  }

  const arma::rowvec ts_num = arma::square(num_sum);
  const arma::rowvec ts_denom = arma::sum(arma::square(denom), 1).t();
  const double ts = n * (ts_num / ts_denom).max();
  return List::create(
    Named("ts") = ts,
    Named("num") = ts_num,
    Named("denom") = ts_denom
  );
}

}  // namespace

//' Fast Zhang--Wang--Shao (2022) q=2 and q=6 statistics
//'
//' This is algebraically equivalent to DMS:::cal_ZWS22(X, q, q).  It removes
//' repeated subvector sums, avoids materializing per-coordinate numerator and
//' denominator cubes, and computes both published q values in one compiled
//' call.  P-values still use the original DMS package calibration samples.
//'
//' @param X Numeric n by p data matrix.
//' @return A list containing q=2 and q=6 statistics plus diagnostic vectors.
// [[Rcpp::export]]
List zws22_fast_stats(const arma::mat& X) {
  const List q2 = zws_one_fast(X, 2, 2);
  const List q6 = zws_one_fast(X, 6, 6);
  return List::create(
    Named("ts2") = q2["ts"],
    Named("ts6") = q6["ts"],
    Named("num2") = q2["num"],
    Named("denom2") = q2["denom"],
    Named("num6") = q6["num"],
    Named("denom6") = q6["denom"]
  );
}
