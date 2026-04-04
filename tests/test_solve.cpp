#include <cassert>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <unistd.h>

#include "SymmetricSparse.h"

// =====================================================================
// Helpers
// =====================================================================

static const char* current_test = nullptr;

void timeout_handler(int) {
    const char prefix[] = "  HANG  ";
    write(STDERR_FILENO, prefix, sizeof(prefix) - 1);
    if (current_test)
        write(STDERR_FILENO, current_test, strlen(current_test));
    write(STDERR_FILENO, "\n", 1);
    _exit(1);
}

void check_close(double a, double b, double tol, const char* msg) {
    if (std::fabs(a - b) > tol) {
        std::cerr << "FAIL " << msg << ": expected " << b
                  << ", got " << a << "\n";
        assert(false);
    }
}

std::vector<int> identity_perm(int n) {
    std::vector<int> p(n);
    for (int i = 0; i < n; ++i) p[i] = i + 1;
    return p;
}

/// Construct, factorize, solve with a given ordering, check solution.
/// order = 0 (identity perm) or -1 (MMD).
void solve_and_check(const char* name, int order,
                     int n,
                     const int* colptr, const int* rowind,
                     const double* nzvals,
                     const double* rhs,
                     const double* x_exact,
                     double tol) {
    current_test = name;
    alarm(5);

    if (order == 0) {
        auto perm = identity_perm(n);
        NgPeytonCpp::SymmetricSparse<double> mat(n, colptr, rowind, 0, &perm[0]);
        mat.ldlTFactorize(colptr, rowind, nzvals);
        std::vector<double> x(n, 0.0);
        mat.solve(rhs, &x[0]);
        alarm(0);
        for (int i = 0; i < n; ++i)
            check_close(x[i], x_exact[i], tol, name);
    } else {
        NgPeytonCpp::SymmetricSparse<double> mat(n, colptr, rowind, -1);
        mat.ldlTFactorize(colptr, rowind, nzvals);
        std::vector<double> x(n, 0.0);
        mat.solve(rhs, &x[0]);
        alarm(0);
        for (int i = 0; i < n; ++i)
            check_close(x[i], x_exact[i], tol, name);
    }
}

/// Run the same test with both identity and MMD orderings.
void run_both(const char* name,
              int n,
              const int* colptr, const int* rowind,
              const double* nzvals,
              const double* rhs,
              const double* x_exact,
              double tol = 1e-12) {
    std::string label_id  = std::string(name) + " (identity)";
    std::string label_mmd = std::string(name) + " (MMD)";

    solve_and_check(label_id.c_str(),  0, n, colptr, rowind, nzvals, rhs, x_exact, tol);
    std::cout << "  PASS  " << label_id << "\n";

    solve_and_check(label_mmd.c_str(), -1, n, colptr, rowind, nzvals, rhs, x_exact, tol);
    std::cout << "  PASS  " << label_mmd << "\n";
}

// =====================================================================
// Test cases
// =====================================================================

// 1x1 scalar: A = [5], rhs = 10, x = 2
void test_1x1() {
    const int n = 1;
    const int colptr[] = {1, 2};
    const int rowind[] = {1};
    const double nzvals[] = {5.0};
    const double rhs[] = {10.0};
    const double x_exact[] = {2.0};
    run_both("1x1", n, colptr, rowind, nzvals, rhs, x_exact);
}

// 3x3 tridiagonal: A = [4 -1 0; -1 4 -1; 0 -1 4], x = {1,1,1}
void test_3x3_tridiagonal() {
    const int n = 3;
    const int colptr[] = {1, 3, 5, 6};
    const int rowind[] = {1, 2, 2, 3, 3};
    const double nzvals[] = {4.0, -1.0, 4.0, -1.0, 4.0};
    const double rhs[] = {3.0, 2.0, 3.0};
    const double x_exact[] = {1.0, 1.0, 1.0};
    run_both("3x3 tridiagonal", n, colptr, rowind, nzvals, rhs, x_exact);
}

// 4x4 diagonal: A = diag(2,3,5,7), x = {2,3,5,7}
void test_4x4_diagonal() {
    const int n = 4;
    const int colptr[] = {1, 2, 3, 4, 5};
    const int rowind[] = {1, 2, 3, 4};
    const double nzvals[] = {2.0, 3.0, 5.0, 7.0};
    const double rhs[] = {4.0, 9.0, 25.0, 49.0};
    const double x_exact[] = {2.0, 3.0, 5.0, 7.0};
    run_both("4x4 diagonal", n, colptr, rowind, nzvals, rhs, x_exact);
}

// 4x4 tridiagonal, x = {1,2,3,4}
void test_4x4_tridiagonal() {
    const int n = 4;
    const int colptr[] = {1, 3, 5, 7, 8};
    const int rowind[] = {1, 2, 2, 3, 3, 4, 4};
    const double nzvals[] = {4.0, -1.0, 4.0, -1.0, 4.0, -1.0, 4.0};
    const double rhs[] = {2.0, 4.0, 6.0, 13.0};
    const double x_exact[] = {1.0, 2.0, 3.0, 4.0};
    run_both("4x4 tridiagonal", n, colptr, rowind, nzvals, rhs, x_exact);
}

// 5x5 arrow: dense first column, x = {1,1,1,1,1}
void test_5x5_arrow() {
    const int n = 5;
    const int colptr[] = {1, 6, 7, 8, 9, 10};
    const int rowind[] = {1, 2, 3, 4, 5, 2, 3, 4, 5};
    const double nzvals[] = {4.0, -1.0, -1.0, -1.0, -1.0,
                             4.0, 4.0, 4.0, 4.0};
    const double rhs[] = {0.0, 3.0, 3.0, 3.0, 3.0};
    const double x_exact[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    run_both("5x5 arrow", n, colptr, rowind, nzvals, rhs, x_exact);
}

// 5x5 dense SPD: all nodes connected, maximum degree, x = {1,1,1,1,1}
void test_5x5_dense() {
    const int n = 5;
    const int colptr[] = {1, 6, 10, 13, 15, 16};
    const int rowind[] = {1, 2, 3, 4, 5,
                          2, 3, 4, 5,
                          3, 4, 5,
                          4, 5,
                          5};
    const double nzvals[] = {10.0, -1.0, -1.0, -1.0, -1.0,
                              10.0, -1.0, -1.0, -1.0,
                              10.0, -1.0, -1.0,
                              10.0, -1.0,
                              10.0};
    const double rhs[] = {6.0, 6.0, 6.0, 6.0, 6.0};
    const double x_exact[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    run_both("5x5 dense", n, colptr, rowind, nzvals, rhs, x_exact);
}

// 2D Laplace 5-point stencil on 8x8 grid (n=64), x = 1 everywhere
void test_laplace2d_8x8() {
    const int nx = 8, ny = 8, n = nx * ny;

    std::vector<int> colptr(n + 1, 0);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            int col = j * nx + i;
            int cnt = 1;
            if (i + 1 < nx) ++cnt;
            if (j + 1 < ny) ++cnt;
            colptr[col + 1] = cnt;
        }
    colptr[0] = 1;
    for (int k = 1; k <= n; ++k) colptr[k] += colptr[k - 1];

    int nnz = colptr[n] - 1;
    std::vector<int> rowind(nnz);
    std::vector<double> nzvals(nnz);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            int col = j * nx + i;
            int pos = colptr[col] - 1;
            rowind[pos] = col + 1; nzvals[pos] = 4.0; ++pos;
            if (i + 1 < nx) { rowind[pos] = col + 2; nzvals[pos] = -1.0; ++pos; }
            if (j + 1 < ny) { rowind[pos] = col + nx + 1; nzvals[pos] = -1.0; ++pos; }
        }

    std::vector<double> rhs(n);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            int nb = 0;
            if (i > 0) ++nb; if (i+1<nx) ++nb; if (j>0) ++nb; if (j+1<ny) ++nb;
            rhs[j * nx + i] = 4.0 - nb;
        }

    std::vector<double> x_exact(n, 1.0);
    run_both("Laplace2D 8x8", n, &colptr[0], &rowind[0], &nzvals[0],
             &rhs[0], &x_exact[0], 1e-10);
}

// =====================================================================
int main() {
    signal(SIGALRM, timeout_handler);
    std::cout << "Running solve tests (identity + MMD)...\n";
    test_1x1();
    test_3x3_tridiagonal();
    test_4x4_diagonal();
    test_4x4_tridiagonal();
    test_5x5_arrow();
    test_5x5_dense();
    test_laplace2d_8x8();
    std::cout << "All tests passed.\n";
    return 0;
}
