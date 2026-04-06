#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "NgPeytonCpp/LDLtSolver.h"
#include "nd2d.h"

using doublecomplex = std::complex<double>;

int main(int argc, char** argv) {
  int i, j, nnodes, nnz;
  int nx = -1, ny = -1, count, node;
  auto order = NgPeytonCpp::Ordering::UserProvided;
  int printa = 0;

  std::vector<int> grid, perm;
  std::vector<int> rowind, colptr;
  std::vector<doublecomplex> nzvals;

  int ia = 1;
  while (ia < argc) {
    if (!strncmp(argv[ia], "-nx", 3)) {
      nx = atoi(&argv[ia][4]);
    } else if (!strncmp(argv[ia], "-ny", 3)) {
      ny = atoi(&argv[ia][4]);
    } else if (!strncmp(argv[ia], "-order", 6)) {
      order = static_cast<NgPeytonCpp::Ordering>(atoi(&argv[ia][7]));
    } else if (!strncmp(argv[ia], "-printa", 7)) {
      printa = atoi(&argv[ia][8]);
      if (printa != 0)
        printa = 1;
    } else {
      std::cerr << " Invalid argument!\n";
      std::cerr << "usage: clap2d5pt -nx=<x> -ny=<y> -order=<n> -printa=1\n";
      return 1;
    }
    ia++;
  }

  if (nx == -1 && ny > 0)
    nx = ny;
  if (ny == -1 && nx > 0)
    ny = nx;

  if (nx < 0)
    nx = 10;
  if (ny < 0)
    ny = 10;

  nnodes = nx * ny;
  grid.resize(nnodes);

  for (i = 0; i < nnodes; i++)
    grid[i] = i;

  // 0-based grid accessor: node at row i, column j
  auto gnode = [&](int i, int j) { return grid[nx * j + i]; };

  // Periodic BC: each node has exactly 4 neighbors
  int nedges = 4 * nnodes;
  doublecomplex dval(10.0, 0.0);
  std::cout << " Periodic boundary condition\n";

  nnz = nedges / 2 + nnodes;

  colptr.resize(nnodes + 1, 0);
  rowind.resize(nnz, 0);
  nzvals.resize(nnz, 0);

  colptr[0] = 0;
  count = 0;
  node = 0;

  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      // Diagonal
      if (printa) {
        std::cout << gnode(i, j) << " " << gnode(i, j) << "  " << real(dval)
                  << " " << imag(dval) << "\n";
      }
      rowind[count] = gnode(i, j);
      nzvals[count] = dval;
      count++;

      // Lower neighbor (i+1, j)
      if (i + 1 < nx) {
        if (printa)
          std::cout << gnode(i + 1, j) << " " << gnode(i, j) << " -1.0 0.0\n";
        rowind[count] = gnode(i + 1, j);
        nzvals[count] = doublecomplex(-1.0, 0.0);
        count++;
      }

      // Periodic wrap: connect i=0 to i=nx-1
      if (i == 0) {
        if (printa)
          std::cout << gnode(nx - 1, j) << " " << gnode(i, j) << " -1.0 0.0\n";
        rowind[count] = gnode(nx - 1, j);
        nzvals[count] = doublecomplex(-1.0, 0.0);
        count++;
      }

      // Right neighbor (i, j+1)
      if (j + 1 < ny) {
        if (printa)
          std::cout << gnode(i, j + 1) << " " << gnode(i, j) << " -1.0 0.0\n";
        rowind[count] = gnode(i, j + 1);
        nzvals[count] = doublecomplex(-1.0, 0.0);
        count++;
      }

      // Periodic wrap: connect j=0 to j=ny-1
      if (j == 0) {
        if (printa)
          std::cout << gnode(i, ny - 1) << " " << gnode(i, j) << " -1.0 0.0\n";
        rowind[count] = gnode(i, ny - 1);
        nzvals[count] = doublecomplex(-1.0, 0.0);
        count++;
      }

      node++;
      colptr[node] = count;
    }
  }
  if (count != nnz) {
    std::cerr << " count = " << count << " nnz = " << nnz << "\n";
    return 1;
  }

  if (order == NgPeytonCpp::Ordering::UserProvided) {
    perm.resize(nnodes);
    nd2d(nx, ny, grid, perm);
  }

  NgPeytonCpp::LDLtSolver<doublecomplex, int> uh(
    nnodes, colptr.data(), rowind.data(), order, perm.data());
  uh.ldlTFactorize(colptr.data(), rowind.data(), nzvals.data());

  // Solve: set x = {0, 1, 2, ...}, compute Ax, solve, check
  std::vector<doublecomplex> x(nnodes), Ax(nnodes, 0.0), y(nnodes, 0.0);
  for (i = 0; i < nnodes; ++i)
    x[i] = doublecomplex(i, 0.0);

  // SpMV: Ax = A*x (lower-triangular CSC, symmetrize)
  for (i = 0; i < nnodes; ++i) {
    for (auto k = colptr[i]; k < colptr[i + 1]; ++k) {
      int r = rowind[k];
      Ax[r] += nzvals[k] * x[i];
      if (r != i)
        Ax[i] += nzvals[k] * x[r];
    }
  }

  uh.solve(Ax.data(), y.data());

  double error = 0.0, norm = 0.0;
  for (i = 0; i < nnodes; ++i) {
    norm += std::abs(x[i]);
    error += std::abs(x[i] - y[i]);
  }
  double relerr = error / norm;
  std::cout << " || x - y ||_1 " << error << "\n";
  std::cout << " || x - y ||_1 / || x ||_1 " << relerr << "\n";

  if (relerr > 1e-12) {
    std::cerr << " FAIL: relative error too large\n";
    return 1;
  }
  return 0;
}
