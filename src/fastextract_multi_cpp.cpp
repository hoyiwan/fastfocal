#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fastextract_multi_cpp(List raster_mats, NumericMatrix coords,
                                    NumericVector x_coords, NumericVector y_coords,
                                    NumericVector res_x, NumericVector res_y) {
  int n_layers = raster_mats.size();
  int npts = coords.nrow();

  NumericMatrix out(npts, n_layers);

  for (int l = 0; l < n_layers; ++l) {
    NumericMatrix mat = raster_mats[l];
    int nrow_raster = mat.nrow();
    int ncol_raster = mat.ncol();

    for (int i = 0; i < npts; ++i) {
      double x = coords(i, 0);
      double y = coords(i, 1);

      int col = (int)((x - x_coords[0]) / res_x[0]);
      int row = (int)((y - y_coords[0]) / res_y[0]);

      if (row >= 0 && row < nrow_raster && col >= 0 && col < ncol_raster) {
        out(i, l) = mat(row, col);
      } else {
        out(i, l) = NA_REAL;
      }
    }
  }

  return out;
}
