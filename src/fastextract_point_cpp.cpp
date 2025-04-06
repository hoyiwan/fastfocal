
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fastextract_point_cpp(NumericMatrix raster_mat, NumericMatrix coords,
                                    NumericVector x_coords, NumericVector y_coords,
                                    NumericVector res_x, NumericVector res_y) {
  int nrow_raster = raster_mat.nrow();
  int ncol_raster = raster_mat.ncol();
  int npts = coords.nrow();
  NumericMatrix out(npts, 1); // for 1-layer raster for now

  for (int i = 0; i < npts; ++i) {
    double x = coords(i, 0);
    double y = coords(i, 1);

    int col = (int)((x - x_coords[0]) / res_x[0]);
    int row = (int)((y - y_coords[0]) / res_y[0]);

    if (row >= 0 && row < nrow_raster && col >= 0 && col < ncol_raster) {
      out(i, 0) = raster_mat(row, col);
    } else {
      out(i, 0) = NA_REAL;
    }
  }

  return out;
}
