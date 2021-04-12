/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <algorithm>
#include <cmath>

namespace mpas {
namespace eqns {

// transpose input 2D matrix
template<typename T>
std::vector<std::vector<T>>
transpose(std::vector<std::vector<T>> & in) {
  std::vector<std::vector<T>> out(in[0].size(), std::vector<T>(in.size()));
  for(std::size_t i = 0; i < in.size(); i++) {
    for(std::size_t j = 0; j < in[0].size(); j++) {
      out[j][i] = in[i][j];
    }
  }
  return out;
}

// matrix multiplication - A x B = C
template<typename T>
std::vector<std::vector<T>>
matMul(std::vector<std::vector<T>> & A, std::vector<std::vector<T>> & B) {
  std::vector<std::vector<T>> C(A.size(), std::vector<T>(B[0].size()));

  flog_assert(A[0].size() == B.size(), "Matrix multiply dimensionality mismatch");

  for(std::size_t row = 0; row < C.size(); row++) {
    for(std::size_t col = 0; col < C[0].size(); col++) {
      for(std::size_t inner = 0; inner < B.size(); inner++) {
        C[row][col] += A[row][inner] * B[inner][col];
      }
    }
  }
  return C;
}

//  subroutine to perform the partial-pivoting gaussian elimination.
//  a(n,n) is the original matrix in the input and transformed matrix
//  plus the pivoting element ratios below the diagonal in the output.
//  indx(n) records the pivoting order.  copyright (c) tao pang 2001.

// translated from src/operators/mpas_matrix_operations.F

template<typename T>
void
elgs(std::vector<std::vector<T>> & a, std::vector<std::size_t> & idx) {
  flog_assert(a.size() == a[0].size(),
    "Gaussian elimination routine designed for square matrices only");
  std::size_t n = a.size();

  // initialize the index
  for(std::size_t i = 0; i < n; i++) {
    idx[i] = i;
  }

  std::vector<T> c(n);
  T c1;

  // find the rescaling factors, one from each row
  for(std::size_t i = 0; i < n; i++) {
    c1 = 0.;
    for(std::size_t j = 0; j < n; j++) {
      c1 = std::max(c1, abs(a[i][j]));
    }
    c[i] = c1;
  }

  // search the pivoting (largest) element from each column
  for(std::size_t j = 0; j < n - 1; j++) {
    std::size_t k;
    T pi1 = 0.0;
    for(std::size_t i = j; i < n; i++) {
      T pi = abs(a[idx[i]][j]) / c[idx[i]];
      if(pi > pi1) {
        pi1 = pi;
        k = i;
      }
    }

    // interchange the rows via indx(n) to record pivoting order
    std::size_t itmp = idx[j];
    idx[j] = idx[k];
    idx[k] = itmp;

    for(std::size_t i = j + 1; i < n; i++) {
      T pj = a[idx[i]][j] / a[idx[j]][j];

      // record pivoting ratios below the diagonal
      a[idx[i]][j] = pj;

      // modify other elements accordingly
      for(std::size_t k = j + 1; k < n; k++) {
        a[idx[i]][k] = a[idx[i]][k] - pj * a[idx[j]][k];
      }
    }
  }
}
//   !!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 4.4   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//   !
//   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//   !                                                                       !
//   ! Please Note:                                                          !
//   !                                                                       !
//   ! (1) This computer program is written by Tao Pang in conjunction with  !
//   !     his book, "An Introduction to Computational Physics," published   !
//   !     by Cambridge University Press in 1997.                            !
//   !                                                                       !
//   ! (2) No warranties, express or implied, are made for this program.     !
//   !                                                                       !
//   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//  subroutine to invert matrix a(n,n) with the inverse stored
//  in x(n,n) in the output.  copyright (c) tao pang 2001.

// translated from src/operators/mpas_matrix_operations.F

template<typename T>
std::vector<std::vector<T>>
migs(std::vector<std::vector<T>> & a) {
  flog_assert(a.size() == a[0].size(),
    "Can only compute inverse for square matrices only");

  std::size_t n = a.size();
  std::vector<std::vector<T>> b(n, std::vector<T>(n));
  std::vector<std::vector<T>> x(n, std::vector<T>(n));

  std::vector<std::size_t> idx(n);

  for(std::size_t j = 0; j < n; j++) {
    for(std::size_t i = 0; i < n; i++) {
      b[i][j] = 0.;
    }
  }
  for(std::size_t i = 0; i < n; i++) {
    b[i][i] = 1.;
  }

  elgs(a, idx);

  for(std::size_t i = 0; i < n - 1; i++) {
    for(std::size_t j = i + 1; j < n; j++) {
      for(std::size_t k = 0; k < n; k++) {
        b[idx[j]][k] = b[idx[j]][k] - a[idx[j]][i] * b[idx[i]][k];
      }
    }
  }

  for(int i = 0; i < n; i++) {
    x[n - 1][i] = b[idx[n - 1]][i] / a[idx[n - 1]][n - 1];
    for(int j = n - 2; j >= 0; j--) {
      x[j][i] = b[idx[j]][i];
      for(int k = j + 1; k < n; k++) {
        x[j][i] = x[j][i] - a[idx[j]][k] * x[k][i];
      }
      x[j][i] = x[j][i] / a[idx[j]][j];
    }
  }
  return x;
}


}}
