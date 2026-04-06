#ifndef LAP2D9PT_H
#define LAP2D9PT_H

#include <cmath>
#include <iostream>
#include <vector>

#include "NgPeytonCpp/LDLtSolver.h"

template <typename Scalar>
int lap2d9pt(int argc, char** argv) {
  int i, j, nnodes, ia;
  int nx = -1, ny = -1, node;

  std::vector<int> grid;
  std::vector<int> rowind, colptr;
  std::vector<Scalar> nzvals;

  ia = 1;
  while (ia < argc) {
    if (!strncmp(argv[ia], "-n", 2)) {
      nx = atoi(&argv[ia][3]);
    } else {
      std::cerr << " Invalid argument!\n";
      std::cerr << "usage: lap2d9pt -n=<x>\n";
      return 1;
    }
    ia++;
  }

  if (nx < 0)
    nx = 11;
  ny = nx;

  nnodes = nx * ny;

  grid.resize(nnodes);
  for (i = 0; i < nnodes; i++)
    grid[i] = i;

  // 0-based grid accessor: node at row i, column j
  auto gnode = [&](int ii, int jj) { return grid[nx * jj + ii]; };

  colptr.resize(nnodes + 1, 0);
  rowind.reserve(9 * nnodes);
  nzvals.reserve(9 * nnodes);

  colptr[0] = 0;
  node = 0;

  Scalar h = Scalar(1) / Scalar(nx + 1);

  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      if (j > 0) {
        if (i > 0) {
          rowind.push_back(gnode(i - 1, j - 1));
          nzvals.push_back(Scalar(-0.25) / (h * h));
        }
        rowind.push_back(gnode(i, j - 1));
        nzvals.push_back(Scalar(-0.5) / (h * h));
        if (i + 1 < nx) {
          rowind.push_back(gnode(i + 1, j - 1));
          nzvals.push_back(Scalar(-0.25) / (h * h));
        }
      }

      if (i > 0) {
        rowind.push_back(gnode(i - 1, j));
        nzvals.push_back(Scalar(-0.5) / (h * h));
      }
      rowind.push_back(gnode(i, j));
      nzvals.push_back(Scalar(3.0) / (h * h));
      if (i + 1 < nx) {
        rowind.push_back(gnode(i + 1, j));
        nzvals.push_back(Scalar(-0.5) / (h * h));
      }

      if (j + 1 < ny) {
        if (i > 0) {
          rowind.push_back(gnode(i - 1, j + 1));
          nzvals.push_back(Scalar(-0.25) / (h * h));
        }
        rowind.push_back(gnode(i, j + 1));
        nzvals.push_back(Scalar(-0.5) / (h * h));
        if (i + 1 < nx) {
          rowind.push_back(gnode(i + 1, j + 1));
          nzvals.push_back(Scalar(-0.25) / (h * h));
        }
      }

      colptr[node + 1] = int(rowind.size());
      node++;
    }
  }

  NgPeytonCpp::LDLtSolver<Scalar, int> LDLt(
    nnodes, colptr.data(), rowind.data(), NgPeytonCpp::Ordering::MMD);
  LDLt.ldlTFactorize(colptr.data(), rowind.data(), nzvals.data());

  std::vector<Scalar> x(nnodes), Ax(nnodes, Scalar(0)), y(nnodes, Scalar(0));
  for (i = 0; i < nnodes; ++i)
    x[i] = Scalar(i);

  for (i = 0; i < nnodes; ++i) {
    Ax[i] = 0;
    for (auto k = colptr[i]; k < colptr[i + 1]; ++k)
      Ax[i] += nzvals[k] * x[rowind[k]];
  }

  LDLt.solve(Ax.data(), y.data());

  double error = 0.0, norm = 0.0;
  for (i = 0; i < nnodes; ++i) {
    norm += double(std::abs(x[i]));
    error += double(std::abs(x[i] - y[i]));
  }
  std::cout << " || x - y ||_1 " << error << "\n";
  std::cout << " || x - y ||_1 / || x ||_1 " << error / norm << "\n";

  return 0;
}

#endif  // LAP2D9PT_H
