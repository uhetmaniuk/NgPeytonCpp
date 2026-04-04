#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "SymmetricSparse.h"

/// Check that |a - b| < tol
void check_close(double a, double b, double tol, const char* msg) {
    if (std::fabs(a - b) > tol) {
        std::cerr << "FAIL " << msg << ": expected " << b
                  << ", got " << a << "\n";
        assert(false);
    }
}

/// Build a 1-based identity permutation of size n
std::vector<int> identity_perm(int n) {
    std::vector<int> p(n);
    for (int i = 0; i < n; ++i) p[i] = i + 1;
    return p;
}

// ---------------------------------------------------------------
// 3x3 tridiagonal SPD matrix (lower-triangular, 1-based CSC)
//
//     A = [  4  -1   0 ]      rhs = [ 3 ]      x = [ 1 ]
//         [ -1   4  -1 ]            [ 2 ]          [ 1 ]
//         [  0  -1   4 ]            [ 3 ]          [ 1 ]
//
// Lower-triangular CSC (1-based):
//   col 1: rows {1,2}  vals {4,-1}
//   col 2: rows {2,3}  vals {4,-1}
//   col 3: rows {3}    vals {4}
// ---------------------------------------------------------------
void test_solve_3x3_tridiagonal() {
    const int n = 3;
    const int colptr[] = {1, 3, 5, 6};
    const int rowind[] = {1, 2, 2, 3, 3};
    const double nzvals[] = {4.0, -1.0, 4.0, -1.0, 4.0};
    auto perm = identity_perm(n);

    NgPeytonCpp::SymmetricSparse<double> mat(n, colptr, rowind, 0, &perm[0]);
    mat.ldlTFactorize(colptr, rowind, nzvals);

    double rhs[] = {3.0, 2.0, 3.0};
    double x[3] = {};
    mat.solve(rhs, x);

    const double tol = 1e-12;
    check_close(x[0], 1.0, tol, "x[0]");
    check_close(x[1], 1.0, tol, "x[1]");
    check_close(x[2], 1.0, tol, "x[2]");
    std::cout << "  PASS  test_solve_3x3_tridiagonal\n";
}

// ---------------------------------------------------------------
// 4x4 diagonal matrix (lower-triangular, 1-based CSC)
//
//     A = diag(2, 3, 5, 7)    rhs = {4, 9, 25, 49}
//     x = {2, 3, 5, 7}
// ---------------------------------------------------------------
void test_solve_4x4_diagonal() {
    const int n = 4;
    const int colptr[] = {1, 2, 3, 4, 5};
    const int rowind[] = {1, 2, 3, 4};
    const double nzvals[] = {2.0, 3.0, 5.0, 7.0};
    auto perm = identity_perm(n);

    NgPeytonCpp::SymmetricSparse<double> mat(n, colptr, rowind, 0, &perm[0]);
    mat.ldlTFactorize(colptr, rowind, nzvals);

    double rhs[] = {4.0, 9.0, 25.0, 49.0};
    double x[4] = {};
    mat.solve(rhs, x);

    const double tol = 1e-12;
    check_close(x[0], 2.0, tol, "x[0]");
    check_close(x[1], 3.0, tol, "x[1]");
    check_close(x[2], 5.0, tol, "x[2]");
    check_close(x[3], 7.0, tol, "x[3]");
    std::cout << "  PASS  test_solve_4x4_diagonal\n";
}

// ---------------------------------------------------------------
// 4x4 tridiagonal SPD matrix (lower-triangular, 1-based CSC)
//
//     A = [  4  -1   0   0 ]
//         [ -1   4  -1   0 ]
//         [  0  -1   4  -1 ]
//         [  0   0  -1   4 ]
//
// A*x with x={1,2,3,4}:
//   row 1:  4*1 - 1*2           =  2
//   row 2: -1*1 + 4*2 - 1*3    =  4
//   row 3:       -1*2 + 4*3 -1*4 =  6
//   row 4:             -1*3 + 4*4 = 13
// ---------------------------------------------------------------
void test_solve_4x4_tridiagonal() {
    const int n = 4;
    const int colptr[] = {1, 3, 5, 7, 8};
    const int rowind[] = {1, 2, 2, 3, 3, 4, 4};
    const double nzvals[] = {4.0, -1.0, 4.0, -1.0, 4.0, -1.0, 4.0};
    auto perm = identity_perm(n);

    NgPeytonCpp::SymmetricSparse<double> mat(n, colptr, rowind, 0, &perm[0]);
    mat.ldlTFactorize(colptr, rowind, nzvals);

    double rhs[] = {2.0, 4.0, 6.0, 13.0};
    double x[4] = {};
    mat.solve(rhs, x);

    const double tol = 1e-12;
    check_close(x[0], 1.0, tol, "x[0]");
    check_close(x[1], 2.0, tol, "x[1]");
    check_close(x[2], 3.0, tol, "x[2]");
    check_close(x[3], 4.0, tol, "x[3]");
    std::cout << "  PASS  test_solve_4x4_tridiagonal\n";
}

// ---------------------------------------------------------------
// 5x5 arrow matrix (lower-triangular, 1-based CSC)
//
//     A = [  4  -1  -1  -1  -1 ]
//         [ -1   4   0   0   0 ]
//         [ -1   0   4   0   0 ]
//         [ -1   0   0   4   0 ]
//         [ -1   0   0   0   4 ]
//
// x = {1,1,1,1,1}, A*x = {0, 3, 3, 3, 3}
// ---------------------------------------------------------------
void test_solve_5x5_arrow() {
    const int n = 5;
    const int colptr[] = {1, 6, 7, 8, 9, 10};
    const int rowind[] = {1, 2, 3, 4, 5, 2, 3, 4, 5};
    const double nzvals[] = {4.0, -1.0, -1.0, -1.0, -1.0,
                             4.0, 4.0, 4.0, 4.0};
    auto perm = identity_perm(n);

    NgPeytonCpp::SymmetricSparse<double> mat(n, colptr, rowind, 0, &perm[0]);
    mat.ldlTFactorize(colptr, rowind, nzvals);

    double rhs[] = {0.0, 3.0, 3.0, 3.0, 3.0};
    double x[5] = {};
    mat.solve(rhs, x);

    const double tol = 1e-12;
    for (int i = 0; i < 5; ++i) {
        check_close(x[i], 1.0, tol, "x[i]");
    }
    std::cout << "  PASS  test_solve_5x5_arrow\n";
}

int main() {
    std::cout << "Running solve tests...\n";
    test_solve_3x3_tridiagonal();
    test_solve_4x4_diagonal();
    test_solve_4x4_tridiagonal();
    test_solve_5x5_arrow();
    std::cout << "All tests passed.\n";
    return 0;
}
