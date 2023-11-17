#ifndef OBJECT_TRACKER_INTERNAL_MATHS_H
#define OBJECT_TRACKER_INTERNAL_MATHS_H

#include <cassert>
#include <iostream>
#include <vector>

namespace tracker {

template <typename T>
void mprint(std::vector<std::vector<T>> A) {
  for (const auto &row : A) {
    for (auto &val : row) {
      std::cout << val << " ";
    }
    std::cout << "\n";
  }
}

template <typename T>
std::vector<std::vector<T>> mmultiply(const std::vector<std::vector<T>> &A,
                                      const std::vector<std::vector<T>> &B) {
  const size_t nA_row = A.size();
  const size_t nA_col = A[0].size();
  const size_t nB_row = B.size();
  const size_t nB_col = B[0].size();

  if (nA_col != nB_row) {
    fprintf(stderr,
            "Number of columns of matrix A must be equal to number of rows of "
            "matrix B\n");
    assert(false);
  }

  std::vector<std::vector<T>> ret(nA_row, std::vector<T>(nB_col, 0));
  for (size_t i = 0; i < nA_row; i++) {
    for (size_t j = 0; j < nB_col; j++) {
      for (size_t k = 0; k < nA_col; k++) {
        ret[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return ret;
}

template <typename T, typename U>
std::vector<std::vector<T>> mmultiply(const std::vector<std::vector<T>> &A,
                                      const U coff) {
  const size_t nA_row = A.size();
  const size_t nA_col = A[0].size();

  std::vector<std::vector<T>> ret(nA_row, std::vector<T>(nA_col, 0));
  for (size_t i = 0; i < nA_row; i++) {
    for (size_t j = 0; j < nA_col; j++) {
      ret[i][j] = A[i][j] * coff;
    }
  }
  return ret;
}

template <typename T, typename... Args>
std::vector<std::vector<T>> mmultiply(const std::vector<std::vector<T>> &A,
                                      const std::vector<std::vector<T>> &B,
                                      Args... args) {
  return mmultiply(A, mmultiply(B, args...));
}

template <typename T>
std::vector<std::vector<T>> madd(const std::vector<std::vector<T>> &A,
                                 const std::vector<std::vector<T>> &B) {
  const size_t nA_row = A.size();
  const size_t nA_col = A[0].size();
  const size_t nB_row = B.size();
  const size_t nB_col = B[0].size();

  if ((nA_row != nB_row) || (nA_col != nB_col)) {
    fprintf(stderr, "Two matrices must have the same size\n");
    assert(false);
  }

  std::vector<std::vector<T>> ret(nA_row, std::vector<T>(nA_col, 0));
  for (size_t i = 0; i < nA_row; i++) {
    for (size_t j = 0; j < nA_col; j++) {
      ret[i][j] = A[i][j] + B[i][j];
    }
  }
  return ret;
}

template <typename T, typename... Args>
std::vector<std::vector<T>> madd(const std::vector<std::vector<T>> &A,
                                 const std::vector<std::vector<T>> &B,
                                 Args... args) {
  return madd(A, madd(B, args...));
}

template <typename T>
std::vector<std::vector<T>> msubtract(const std::vector<std::vector<T>> &A,
                                      const std::vector<std::vector<T>> &B) {
  return madd(A, mmultiply(B, -1.0f));
}

template <typename T>
std::vector<std::vector<T>> mtranspose(const std::vector<std::vector<T>> &A) {
  const size_t nA_row = A.size();
  const size_t nA_col = A[0].size();

  std::vector<std::vector<T>> ret(nA_col, std::vector<T>(nA_row, 0));
  for (size_t i = 0; i < nA_row; i++) {
    for (size_t j = 0; j < nA_col; j++) {
      ret[j][i] = A[i][j];
    }
  }
  return ret;
}

template <typename T = float>
std::vector<std::vector<T>> midentity(size_t size) {
  std::vector<std::vector<T>> ret(size, std::vector<T>(size, 0));
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      if (i == j) {
        ret[i][j] = static_cast<T>(1);
      }
    }
  }
  return ret;
}

}  // namespace tracker

#endif  // OBJECT_TRACKER_INTERNAL_MATHS_H
