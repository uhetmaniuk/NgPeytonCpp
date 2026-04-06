#include <cassert>
#include <iostream>

// NOTE: This test exercises the internal isFullRepresentation() function
// directly. It uses 1-based CSC indices (the internal convention), not the
// 0-based indices of the public API.
#include "details/Utilities.h"

using NgPeytonCpp::details::isFullRepresentation;

// ---------------------------------------------------------------
// Helper: 3x3 lower-triangular only (1-based CSC)
//
//     [ 1  .  . ]
//     [ 2  4  . ]
//     [ 3  5  6 ]
//
// Lower-triangular CSC (1-based):
//   col 1: rows {1,2,3}   col 2: rows {2,3}   col 3: rows {3}
//   colptr = {1, 4, 6, 7}
//   rowind = {1, 2, 3,  2, 3,  3}
// ---------------------------------------------------------------
void test_lower_triangular_only() {
  const int n = 3;
  const int colptr[] = {1, 4, 6, 7};
  const int rowind[] = {1, 2, 3, 2, 3, 3};

  bool result = isFullRepresentation(n, colptr, rowind);
  assert(!result && "Lower-triangular matrix must NOT be detected as full");
  std::cout << "  PASS  test_lower_triangular_only\n";
}

// ---------------------------------------------------------------
// Helper: 3x3 full representation (1-based CSC)
//
//     [ 1  2  3 ]
//     [ 2  4  5 ]
//     [ 3  5  6 ]
//
// Full CSC (1-based):
//   col 1: rows {1,2,3}   col 2: rows {1,2,3}   col 3: rows {1,2,3}
//   colptr = {1, 4, 7, 10}
//   rowind = {1, 2, 3,  1, 2, 3,  1, 2, 3}
// ---------------------------------------------------------------
void test_full_representation() {
  const int n = 3;
  const int colptr[] = {1, 4, 7, 10};
  const int rowind[] = {1, 2, 3, 1, 2, 3, 1, 2, 3};

  bool result = isFullRepresentation(n, colptr, rowind);
  assert(result && "Full-representation matrix must be detected as full");
  std::cout << "  PASS  test_full_representation\n";
}

// ---------------------------------------------------------------
// Diagonal-only matrix (n=3): every column has only its diagonal
//   colptr = {1, 2, 3, 4}
//   rowind = {1, 2, 3}
//
// No off-diagonal entries, so the function should return false.
// ---------------------------------------------------------------
void test_diagonal_only() {
  const int n = 3;
  const int colptr[] = {1, 2, 3, 4};
  const int rowind[] = {1, 2, 3};

  bool result = isFullRepresentation(n, colptr, rowind);
  assert(!result && "Diagonal-only matrix must NOT be detected as full");
  std::cout << "  PASS  test_diagonal_only\n";
}

// ---------------------------------------------------------------
// Single element matrix (n=1): trivially lower-triangular
//   colptr = {1, 2}
//   rowind = {1}
// ---------------------------------------------------------------
void test_single_element() {
  const int n = 1;
  const int colptr[] = {1, 2};
  const int rowind[] = {1};

  bool result = isFullRepresentation(n, colptr, rowind);
  assert(!result && "1x1 matrix must NOT be detected as full");
  std::cout << "  PASS  test_single_element\n";
}

// ---------------------------------------------------------------
// 9x9 matrix: 6 diagonal columns, then lower-triangular (cols 7-9)
//
//   col 1: {1}  col 2: {2}  ...  col 6: {6}
//   col 7: {7,8,9}   col 8: {8,9}   col 9: {9}
//
//   colptr = {1, 2, 3, 4, 5, 6, 7, 10, 12, 13}
//   rowind = {1, 2, 3, 4, 5, 6,  7, 8, 9,  8, 9,  9}
// ---------------------------------------------------------------
void test_diagonal_then_lower_triangular() {
  const int n = 9;
  const int colptr[] = {1, 2, 3, 4, 5, 6, 7, 10, 12, 13};
  const int rowind[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 9, 9};

  bool result = isFullRepresentation(n, colptr, rowind);
  assert(!result && "Diagonal-then-lower-tri must NOT be detected as full");
  std::cout << "  PASS  test_diagonal_then_lower_triangular\n";
}

// ---------------------------------------------------------------
// 9x9 matrix: 6 diagonal columns, then full-representation (cols 7-9)
//
//   col 1: {1}  col 2: {2}  ...  col 6: {6}
//   col 7: {7,8,9}   col 8: {7,8,9}   col 9: {7,8,9}
//
//   colptr = {1, 2, 3, 4, 5, 6, 7, 10, 13, 16}
//   rowind = {1, 2, 3, 4, 5, 6,  7, 8, 9,  7, 8, 9,  7, 8, 9}
// ---------------------------------------------------------------
void test_diagonal_then_full_representation() {
  const int n = 9;
  const int colptr[] = {1, 2, 3, 4, 5, 6, 7, 10, 13, 16};
  const int rowind[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 7, 8, 9, 7, 8, 9};

  bool result = isFullRepresentation(n, colptr, rowind);
  assert(result && "Diagonal-then-full must be detected as full");
  std::cout << "  PASS  test_diagonal_then_full_representation\n";
}

// ---------------------------------------------------------------
// 9x9 matrix: lower-triangular (cols 1-3), then 6 diagonal columns
//
//   col 1: {1,2,3}   col 2: {2,3}   col 3: {3}
//   col 4: {4}  col 5: {5}  ...  col 9: {9}
//
//   colptr = {1, 4, 6, 7, 8, 9, 10, 11, 12, 13}
//   rowind = {1, 2, 3,  2, 3,  3,  4, 5, 6, 7, 8, 9}
// ---------------------------------------------------------------
void test_lower_triangular_then_diagonal() {
  const int n = 9;
  const int colptr[] = {1, 4, 6, 7, 8, 9, 10, 11, 12, 13};
  const int rowind[] = {1, 2, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9};

  bool result = isFullRepresentation(n, colptr, rowind);
  assert(!result && "Lower-tri-then-diagonal must NOT be detected as full");
  std::cout << "  PASS  test_lower_triangular_then_diagonal\n";
}

// ---------------------------------------------------------------
// 9x9 matrix: full-representation (cols 1-3), then 6 diagonal columns
//
//   col 1: {1,2,3}   col 2: {1,2,3}   col 3: {1,2,3}
//   col 4: {4}  col 5: {5}  ...  col 9: {9}
//
//   colptr = {1, 4, 7, 10, 11, 12, 13, 14, 15, 16}
//   rowind = {1, 2, 3,  1, 2, 3,  1, 2, 3,  4, 5, 6, 7, 8, 9}
// ---------------------------------------------------------------
void test_full_representation_then_diagonal() {
  const int n = 9;
  const int colptr[] = {1, 4, 7, 10, 11, 12, 13, 14, 15, 16};
  const int rowind[] = {1, 2, 3, 1, 2, 3, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  bool result = isFullRepresentation(n, colptr, rowind);
  assert(result && "Full-then-diagonal must be detected as full");
  std::cout << "  PASS  test_full_representation_then_diagonal\n";
}

int main() {
  std::cout << "Running isFullRepresentation tests...\n";
  test_lower_triangular_only();
  test_full_representation();
  test_diagonal_only();
  test_single_element();
  test_diagonal_then_lower_triangular();
  test_diagonal_then_full_representation();
  test_lower_triangular_then_diagonal();
  test_full_representation_then_diagonal();
  std::cout << "All tests passed.\n";
  return 0;
}
