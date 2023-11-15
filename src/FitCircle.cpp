#define STRICT_R_HEADERS
#include <Rcpp.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

struct Point {
  double x;
  double y;
};

// Function to fit a circle to a set of data points
void fitCircle(const std::vector<Point>& data, double& centerX, double& centerY, double& radius) {
  int n = data.size();

  // Create a design matrix
  Eigen::MatrixXd A(n, 3);
  Eigen::VectorXd b(n);

  for (int i = 0; i < n; i++) {
    A(i, 0) = 2 * data[i].x;
    A(i, 1) = 2 * data[i].y;
    A(i, 2) = -1;
    b(i) = data[i].x * data[i].x + data[i].y * data[i].y;
  }

  // Solve for the circle parameters
  Eigen::Vector3d x = A.fullPivLu().solve(b);

  centerX = x(0);
  centerY = x(1);
  radius = std::sqrt(x(0) * x(0) + x(1) * x(1) - x(2));
}

int main() {
  // Sample data (replace with your actual data)
  std::vector<Point> data = {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}};

  double centerX, centerY, radius;
  fitCircle(data, centerX, centerY, radius);

  // Print the results
  std::cout << "Circle Center: (" << centerX << ", " << centerY << ")\n";
  std::cout << "Circle Radius: " << radius << std::endl;

  return 0;
}
