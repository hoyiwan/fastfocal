
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fastfocal_cpp(NumericMatrix mat,
                            double res_x, double res_y,
                            double radius,
                            std::string stat = "mean") {
  int nrows = mat.nrow();
  int ncols = mat.ncol();
  NumericMatrix out(nrows, ncols);
  double radius_sq = radius * radius;

  for (int row_c = 0; row_c < nrows; ++row_c) {
    for (int col_c = 0; col_c < ncols; ++col_c) {
      std::vector<double> vals;

      for (int dr = -std::ceil(radius / res_y); dr <= std::ceil(radius / res_y); ++dr) {
        for (int dc = -std::ceil(radius / res_x); dc <= std::ceil(radius / res_x); ++dc) {
          int row = row_c + dr;
          int col = col_c + dc;

          double dx = dc * res_x;
          double dy = dr * res_y;
          double dist_sq = dx * dx + dy * dy;

          if (dist_sq <= radius_sq &&
              row >= 0 && row < nrows && col >= 0 && col < ncols) {
            double val = mat(row, col);
            if (!NumericVector::is_na(val)) {
              vals.push_back(val);
            }
          }
        }
      }

      if (vals.empty()) {
        out(row_c, col_c) = NA_REAL;
        continue;
      }

      if (stat == "mean") {
        double sum = std::accumulate(vals.begin(), vals.end(), 0.0);
        out(row_c, col_c) = sum / vals.size();
      } else if (stat == "min") {
        out(row_c, col_c) = *std::min_element(vals.begin(), vals.end());
      } else if (stat == "max") {
        out(row_c, col_c) = *std::max_element(vals.begin(), vals.end());
      } else if (stat == "sum") {
        out(row_c, col_c) = std::accumulate(vals.begin(), vals.end(), 0.0);
      } else if (stat == "sd") {
        double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
        double var = 0.0;
        for (double v : vals) var += (v - mean) * (v - mean);
        out(row_c, col_c) = std::sqrt(var / (vals.size() - 1));
      } else if (stat == "median") {
        std::nth_element(vals.begin(), vals.begin() + vals.size() / 2, vals.end());
        double med = vals[vals.size() / 2];
        if (vals.size() % 2 == 0) {
          auto max_low = *std::max_element(vals.begin(), vals.begin() + vals.size() / 2);
          med = (med + max_low) / 2.0;
        }
        out(row_c, col_c) = med;
      } else if (stat == "range") {
        double min = *std::min_element(vals.begin(), vals.end());
        double max = *std::max_element(vals.begin(), vals.end());
        out(row_c, col_c) = max - min;
      } else if (stat == "p25" || stat == "p75") {
        std::sort(vals.begin(), vals.end());
        double p = (stat == "p25") ? 0.25 : 0.75;
        size_t idx = std::floor(p * (vals.size() - 1));
        out(row_c, col_c) = vals[idx];
      } else {
        out(row_c, col_c) = NA_REAL;
      }
    }
  }

  return out;
}
