/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 Inria
 * Author: Emmanuel Benazera <emmanuel.benazera@lri.fr>
 *
 * This file is part of libcmaes.
 *
 * libcmaes is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libcmaes is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libcmaes.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef EO_MATRIX_H
#define EO_MATRIX_H

#include <algorithm>
#include <Eigen/Dense>
#include <stdlib.h>

typedef Eigen::MatrixXd dMat; // declares a column-major non-sparse matrix type of double
typedef Eigen::VectorXd dVec; // declares a vector of double.

#include <unsupported/Eigen/MatrixFunctions>

inline void removeRow(dMat& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();

  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

  matrix.conservativeResize(numRows,numCols);
}

inline void removeColumn(dMat& matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;

  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

  matrix.conservativeResize(numRows,numCols);
}

inline void removeElement(dVec &vec, unsigned int k)
{
  if (k >= vec.size())
    return;
  std::copy(vec.data()+k+1,vec.data()+vec.size(),vec.data()+k);
  vec.conservativeResize(vec.size()-1);
}

inline void addElement(dVec &vec, unsigned int k, const double &xk)
{
  if (k >= vec.size()+1)
    return;
  vec.conservativeResize(vec.size()+1);
  std::copy(vec.data()+k,vec.data()+vec.size()-1,vec.data()+k+1);
  vec[k] = xk;
}

#endif
