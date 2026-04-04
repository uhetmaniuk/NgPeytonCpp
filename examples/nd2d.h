#ifndef ND2D_H
#define ND2D_H

#include <cmath>
#include <vector>

#define MESHMIN 3
#define mesh_(i, j) mesh[nx * (j) + (i)]

int nd2d(int nx, int ny, std::vector<int>& mesh, std::vector<int>& perm) {
  int status = 0;
  int i, j, count;
  int nxtop, nxbot, nyleft, nyright;
  int midrow, midcol;

  if (nx <= MESHMIN && ny <= MESHMIN) {
    /* at a leaf, walk thru in column major */
    count = 0;
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++) {
        perm[count] = mesh_(i, j);
        count++;
      }
  } else if (nx > ny) {
    /* cut mesh in the middle row */
    midrow = (int)ceil((double)nx / 2.0) - 1;  // 0-based index of separator

    /* order the top half (rows 0..midrow-1) */
    nxtop = midrow;
    if (nxtop > 0) {
      std::vector<int> meshtop(nxtop * ny);
      std::vector<int> ptop(nxtop * ny, 0);
      for (j = 0; j < ny; j++)
        for (i = 0; i < nxtop; i++)
          meshtop[nxtop * j + i] = mesh_(i, j);

      nd2d(nxtop, ny, meshtop, ptop);
      for (i = 0; i < nxtop * ny; i++)
        perm[i] = ptop[i];
    }

    /* order the bottom half (rows midrow+1..nx-1) */
    nxbot = nx - nxtop - 1;
    if (nxbot > 0) {
      std::vector<int> meshbot(nxbot * ny);
      std::vector<int> pbot(nxbot * ny, 0);
      for (j = 0; j < ny; j++)
        for (i = 0; i < nxbot; i++)
          meshbot[nxbot * j + i] = mesh_(midrow + 1 + i, j);

      nd2d(nxbot, ny, meshbot, pbot);
      for (i = 0; i < nxbot * ny; i++)
        perm[nxtop * ny + i] = pbot[i];
    }

    /* append the separator (row midrow) */
    for (j = 0; j < ny; j++)
      perm[(nxtop + nxbot) * ny + j] = mesh_(midrow, j);

  } else {
    /* cut mesh in the middle column */
    midcol = (int)ceil((double)ny / 2.0) - 1;  // 0-based index of separator

    /* order the left half (cols 0..midcol-1) */
    nyleft = midcol;
    if (nyleft > 0) {
      std::vector<int> meshleft(nx * nyleft);
      std::vector<int> pleft(nx * nyleft, 0);
      for (j = 0; j < nyleft; j++)
        for (i = 0; i < nx; i++)
          meshleft[j * nx + i] = mesh_(i, j);

      nd2d(nx, nyleft, meshleft, pleft);
      for (i = 0; i < nyleft * nx; i++)
        perm[i] = pleft[i];
    }

    /* order the right half (cols midcol+1..ny-1) */
    nyright = ny - nyleft - 1;
    if (nyright > 0) {
      std::vector<int> meshright(nx * nyright);
      std::vector<int> pright(nx * nyright, 0);
      for (j = 0; j < nyright; j++)
        for (i = 0; i < nx; i++)
          meshright[j * nx + i] = mesh_(i, midcol + 1 + j);

      nd2d(nx, nyright, meshright, pright);
      for (i = 0; i < nyright * nx; i++)
        perm[nyleft * nx + i] = pright[i];
    }

    /* append the separator (col midcol) */
    for (i = 0; i < nx; i++)
      perm[(nyleft + nyright) * nx + i] = mesh_(i, midcol);
  }
  return status;
}

#undef mesh_
#undef MESHMIN

#endif  // ND2D_H
