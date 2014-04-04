
#ifndef IO_MATRIX_H
#define IO_MATRIX_H

#include <Eigen/Sparse>
#include <stdlib.h>

typedef Eigen::MatrixXd dMat; // declares a column-major non-sparse matrix type of double
typedef Eigen::VectorXd dVec; // declares a vector of double.

using namespace Eigen;

#include <unsupported/Eigen/MatrixFunctions>

#endif
