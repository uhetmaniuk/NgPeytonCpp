#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "SymmetricSparse.h"
#include "nd2d.h"

using doublecomplex = std::complex<double>;

#define mesh(i, j) mesh[nx * ((j)-1) + (i)-1]

int main(int argc, char **argv) {

  int i, j, nnodes, nedges, ia, nnz;
  int nx = -1, ny = -1, count, node;
  int token, order = 0;
  int Lnnz;
  double t0 = 0.0, t1 = 0.0, errmax = 0.0;
  doublecomplex dval;
  int pbc = 1, chkerr = 1;
  int printa = 0;
  int dumpL = 0;

  std::vector<int> mesh, perm;
  std::vector<int> rowind, colptr;
  std::vector<int> rowind_inva, colptr_inva;
  std::vector<doublecomplex> nzvals, diag, diag2, inva;

  ia = 1;
  while (ia < argc) {
    if (!strncmp(argv[ia], "-nx", 3)) {
      nx = atoi(&argv[ia][4]);
    } else if (!strncmp(argv[ia], "-ny", 3)) {
      ny = atoi(&argv[ia][4]);
    } else if (!strncmp(argv[ia], "-order", 6)) {
      order = atoi(&argv[ia][7]);
    } else if (!strncmp(argv[ia], "-printa", 7)) {
      printa = atoi(&argv[ia][8]);
      if (printa != 0)
        printa = 1;
    } else {
      std::cerr << " Invalid argument!\n";
      std::cerr << "usage: clap2d5pt -nx=<x> -ny=<y> -order=<n> "
                   "-chkerr=1 -printa=1\n";
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

  nx = 5;
  ny = 7;
  printa = 1;
  chkerr = 1;

  nnodes = nx * ny;
  mesh.resize(nnodes);

  for (i = 0; i < nnodes; i++)
    mesh[i] = i + 1;

  /* first pass to count the number of edges */
  if (pbc) {
    /* Periodic BC */
    nedges = 4 * nnodes;
  } else {
    /* Dirichlet BC */
    nedges = 0;
    for (j = 1; j <= ny; j++) {
      for (i = 1; i <= nx; i++) {
        if (j < ny)
          nedges++;
        if (j > 1)
          nedges++;
        if (i < nx)
          nedges++;
        if (i > 1)
          nedges++;
      }
    }
  }
  /* print the matrix dimension and number of nonzeros */
  nnz = nedges / 2 + nnodes;

  colptr.resize(nnodes + 1, 0);
  rowind.resize(nnz, 0);
  nzvals.resize(nnz, 0);

  colptr[0] = 1;
  count = 0;
  node = 0;

  if (pbc) {
    dval = std::complex<double>(10.0, 0.0);
    std::cout << " Periodic boundary condition\n";
  } else {
    dval = std::complex<double>(4.0, 0.0);
    std::cout << " Dirichlet boundary condition\n";
  }

  for (j = 1; j <= ny; j++) {
    for (i = 1; i <= nx; i++) {
      /* diagonal */
      if (printa) {
        printf("%d %d  %8.2e %8.2e\n", mesh(i, j), mesh(i, j), real(dval),
               imag(dval));
      }

      rowind[count] = mesh(i, j);
      nzvals[count] = dval;
      count++;

      /* lower */
      if (i < nx) {
        if (printa)
          printf("%d %d -1.0 0.0\n", mesh(i + 1, j), mesh(i, j));

        rowind[count] = mesh(i + 1, j);
        nzvals[count] = std::complex<double>(-1.0, 0.0);
        count++;
      }

      if (pbc) {
        /* bottom of the mesh */
        if (i == 1) {
          if (printa)
            printf("%d %d -1.0 0.0\n", mesh(nx, j), mesh(i, j));

          rowind[count] = mesh(nx, j);
          nzvals[count] = std::complex<double>(-1.0, 0.0);
          count++;
        }
      }

      /* right */
      if (j < ny) {
        if (printa)
          printf("%d %d -1.0 0.0\n", mesh(i, j + 1), mesh(i, j));

        rowind[count] = mesh(i, j + 1);
        nzvals[count] = std::complex<double>(-1.0, 0.0);
        count++;
      }

      if (pbc) {
        /* right end of the mesh */
        if (j == 1) {
          if (printa)
            printf("%d %d -1.0 0.0\n", mesh(i, ny), mesh(i, j));

          rowind[count] = mesh(i, ny);
          nzvals[count] = std::complex<double>(-1.0, 0.0);
          count++;
        }
      }
      node++;
      colptr[node] = count + 1;
    }
  }
  if (count != nnz) {
    std::cerr << " count = " << count << " nnz = " << nnz << "\n";
    return 1;
  }

  token = 0;

  if (order == 0) {
    perm.resize(nnodes);
    nd2d(nx, ny, mesh, perm);
  }

  std::vector<int> u_cptr(colptr), u_rowind(rowind);
  std::vector<doublecomplex> u_val(nzvals);
  std::vector<int> u_perm(perm);

  NgPeytonCpp::SymmetricSparse<doublecomplex> uh(nnodes, u_cptr.data(),
                                                 u_rowind.data(), order, u_perm.data());
  uh.ldlTFactorize(&u_cptr[0], &u_rowind[0], &u_val[0]);

}
