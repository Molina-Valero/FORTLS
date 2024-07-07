#define STRICT_R_HEADERS
#include <Rcpp.h>
using namespace Rcpp;

// Function to fit a circle through three points
// Input: a 3x2 matrix representing three points (x, y)
// Output: a list containing center (x, y) and radius
// Adapted from your original R function .fit_circle
// Note: This is a simple implementation and may need further optimizations
//       depending on your specific requirements.

// [[Rcpp::export]]
List fit_circle_cpp(NumericMatrix points) {
  double x1 = points(0, 0);
  double y1 = points(0, 1);
  double x2 = points(1, 0);
  double y2 = points(1, 1);
  double x3 = points(2, 0);
  double y3 = points(2, 1);

  double D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
  double Ux = ((x1 * x1 + y1 * y1) * (y2 - y3) + (x2 * x2 + y2 * y2) * (y3 - y1) + (x3 * x3 + y3 * y3) * (y1 - y2)) / D;
  double Uy = ((x1 * x1 + y1 * y1) * (x3 - x2) + (x2 * x2 + y2 * y2) * (x1 - x3) + (x3 * x3 + y3 * y3) * (x2 - x1)) / D;

  double radius = sqrt((x1 - Ux) * (x1 - Ux) + (y1 - Uy) * (y1 - Uy));

  return List::create(_["center"] = NumericVector::create(Ux, Uy),
                      _["radius"] = radius);
}
