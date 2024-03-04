#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "SymmetricSparse.h"
#include "nd2d.h"

using Scalar = double;

#define mesh(i, j) mesh[nx * (j) + (i)]

int main(int argc, char **argv) {

    int i, j, nnodes, nedges, ia, nnz;
    int nx = -1, ny = -1, count, node;
    int order = 1;

    std::vector<int> mesh, perm;
    std::vector<int> rowind, colptr;
    std::vector<Scalar> nzvals;

    ia = 1;
    while (ia < argc) {
        if (!strncmp(argv[ia], "-nx", 3)) {
            nx = atoi(&argv[ia][4]);
        } else if (!strncmp(argv[ia], "-ny", 3)) {
            ny = atoi(&argv[ia][4]);
        } else if (!strncmp(argv[ia], "-order", 6)) {
            order = atoi(&argv[ia][7]);
        } else {
            std::cerr << " Invalid argument!\n";
            std::cerr << "usage: clap2d5pt -nx=<x> -ny=<y> -order=<n> \n";
            return 1;
        }
        ia++;
    }

    if (nx == -1 && ny > 0)
        nx = ny;
    if (ny == -1 && nx > 0)
        ny = nx;

    if (nx < 0)
        nx = 11;
    if (ny < 0)
        ny = 9;

    nnodes = nx * ny;

    mesh.resize(nnodes);
    for (i = 0; i < nnodes; i++) {
        mesh[i] = i;
    }

    /* first pass to count the number of edges */
    /* Dirichlet BC */
    nedges = 0;
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            if (j + 1 < ny)
                nedges++;
            if (j > 0)
                nedges++;
            if (i + 1 < nx)
                nedges++;
            if (i > 0)
                nedges++;
        }
    }
    /* print the matrix dimension and number of nonzeros */
    nnz = nedges / 2 + nnodes;

    colptr.resize(nnodes + 1, 0);
    rowind.resize(nnz, 0);
    nzvals.resize(nnz, 0);

    colptr[0] = 0;
    count = 0;
    node = 0;

    Scalar hx = Scalar(1) / Scalar(nx + 1);
    Scalar hy = Scalar(1) / Scalar(ny + 1);

    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            /* diagonal */
            rowind[count] = mesh(i, j);
            nzvals[count] = Scalar(2) / (hx * hx) + Scalar(2) / (hy * hy);
            count++;

            /* lower */
            if (i + 1 < nx) {
                rowind[count] = mesh(i + 1, j);
                nzvals[count] = Scalar(-1.0) / (hx * hx);
                count++;
            }

            /* right */
            if (j + 1 < ny) {
                rowind[count] = mesh(i, j + 1);
                nzvals[count] = Scalar(-1.0) / (hy * hy);
                count++;
            }

            colptr[node + 1] = count;
            node++;
        }
    }

    if (count != nnz) {
        std::cerr << " count = " << count << " nnz = " << nnz << "\n";
        return 2;
    }

    if (order == 0) {
        perm.resize(nnodes);
        nd2d(nx, ny, mesh, perm);
    }

    NgPeytonCpp::SymmetricSparse<Scalar> LDLt(nnodes, colptr.data(),
                                              rowind.data(), nzvals.data(), order, perm.data());

    std::vector<Scalar> x(nnodes), Ax(nnodes, Scalar(0)), y(nnodes, Scalar(0));
    for (i = 0; i < nnodes; ++i) {
        x[i] = Scalar(std::rand()) / Scalar(RAND_MAX);
    }

    for (i = 0; i < nnodes; ++i) {
        for (auto k = colptr[i]; k < colptr[i + 1]; ++k) {
            Ax[i] += nzvals[k] * x[rowind[k]];
            if (rowind[k] != i) {
                Ax[rowind[k]] += nzvals[k] * x[i];
            }
        }
    }

    LDLt.solve(Ax.data(), y.data());
    double error = 0.0, norm = 0.0;
    for (i = 0; i < nnodes; ++i) {
        norm += std::abs(x[i]);
        error += std::abs(x[i] - y[i]);
    }
    std::cout << " || x - y || " << error << std::endl;
    std::cout << " || x - y || / || x || " << error / norm << "\n";

    return 0;
}
