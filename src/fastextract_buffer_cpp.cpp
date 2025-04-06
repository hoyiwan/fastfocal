
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fastextract_buffer_cpp(NumericMatrix mat, NumericMatrix coords,
                                     double x_min, double y_min,
                                     double res_x, double res_y,
                                     double radius, std::string stat,
                                     bool na_rm = true, std::string window = "circular") {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  int npts = coords.nrow();
  NumericMatrix out(npts, 1);
  
  double radius_sq = radius * radius;
  double sigma = radius / 2.0;  // for Gaussian falloff
  double two_sigma_sq = 2.0 * sigma * sigma;
  
  for (int i = 0; i < npts; ++i) {
    double x = coords(i, 0);
    double y = coords(i, 1);
    
    int center_col = std::floor((x - x_min) / res_x);
    int center_row = nrow - 1 - std::floor((y - y_min) / res_y);
    
    std::vector<double> values;
    std::vector<double> weights;
    
    for (int dx = -std::ceil(radius); dx <= std::ceil(radius); ++dx) {
      for (int dy = -std::ceil(radius); dy <= std::ceil(radius); ++dy) {
        int row = center_row + dy;
        int col = center_col + dx;
        
        if (row >= 0 && row < nrow && col >= 0 && col < ncol) {
          bool include = false;
          double weight = 1.0;
          
          double dist_sq = dx * dx + dy * dy;
          
          if (window == "circular") {
            include = dist_sq <= radius_sq;
          } else if (window == "rectangular") {
            include = (std::abs(dx) <= radius && std::abs(dy) <= radius);
          } else if (window == "gaussian") {
            include = dist_sq <= radius_sq;
            weight = std::exp(-dist_sq / two_sigma_sq);
          }
          
          if (include) {
            double val = mat(row, col);
            
            if (!Rcpp::NumericVector::is_na(val)) {
              values.push_back(val);
              weights.push_back(weight);
            } else if (!na_rm) {
              values.clear();
              values.push_back(NA_REAL);
              goto finish;
            }
          }
        }
      }
    }
    
    finish:
      if (values.empty()) {
        out(i, 0) = NA_REAL;
      } else {
        if (stat == "mean") {
          if (window == "gaussian") {
            double weighted_sum = 0.0, weight_total = 0.0;
            for (size_t j = 0; j < values.size(); ++j) {
              weighted_sum += values[j] * weights[j];
              weight_total += weights[j];
            }
            out(i, 0) = weighted_sum / weight_total;
          } else {
            out(i, 0) = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
          }
        } else if (stat == "sum") {
          out(i, 0) = std::accumulate(values.begin(), values.end(), 0.0);
        } else if (stat == "min") {
          out(i, 0) = *std::min_element(values.begin(), values.end());
        } else if (stat == "max") {
          out(i, 0) = *std::max_element(values.begin(), values.end());
        } else if (stat == "sd") {
          double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
          double sq_sum = 0.0;
          for (double v : values) sq_sum += (v - mean) * (v - mean);
          out(i, 0) = std::sqrt(sq_sum / (values.size() - 1));
        } else if (stat == "median" || stat == "p50") {
          std::sort(values.begin(), values.end());
          size_t mid = values.size() / 2;
          out(i, 0) = (values.size() % 2 == 0) ? 
          (values[mid - 1] + values[mid]) / 2.0 : 
            values[mid];
        } else if (stat == "range") {
          out(i, 0) = *std::max_element(values.begin(), values.end()) - *std::min_element(values.begin(), values.end());
        } else if (stat == "p25") {
          std::sort(values.begin(), values.end());
          size_t idx = static_cast<size_t>(0.25 * values.size());
          out(i, 0) = values[idx];
        } else if (stat == "p75") {
          std::sort(values.begin(), values.end());
          size_t idx = static_cast<size_t>(0.75 * values.size());
          out(i, 0) = values[idx];
        } else {
          out(i, 0) = NA_REAL; // Unknown stat
        }
      }
  }
  
  return out;
}
