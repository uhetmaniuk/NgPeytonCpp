#ifndef NGPEYTONCPP_DETAILS_BLAS_H
#define NGPEYTONCPP_DETAILS_BLAS_H

#include "NgPeytonCpp/LDLtSolver.h"

namespace NgPeytonCpp { namespace details { namespace f2c {

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     DSCAL .... SCALE A VECTOR                     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE COMPUTES A <-- AX, WHERE A IS A */
/*               SCALAR AND X IS A VECTOR. */

/*     INPUT PARAMETERS - */
/*        N - LENGTH OF THE VECTOR X. */
/*        A - SCALAR MULIPLIER. */
/*        X - VECTOR TO BE SCALED. */

/*     OUTPUT PARAMETERS - */
/*        X - REPLACED BY THE SCALED VECTOR, AX. */

/* *********************************************************************** */

template <typename Scalar>
void scal(int64_t n, Scalar a, Scalar* x) {
  for (int64_t i = 0; i < n; ++i) {
    x[i] = a * x[i];
  }
} /* scal */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     SMXPY1 .... MATRIX-VECTOR MULTIPLY            ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS A MATRIX-VECTOR MULTIPLY, */
/*               Y = Y + AX, ASSUMING DATA STRUCTURES USED IN */
/*               RECENTLY DEVELOPED SPARSE CHOLESKY CODES.  THE */
/*               '1' SIGNIFIES NO LOOP UNROLLING, I.E., */
/*               LOOP-UNROLLING TO LEVEL 1. */

/*     INPUT PARAMETERS - */
/*        M      - NUMBER OF ROWS. */
/*        N      - NUMBER OF COLUMNS. */
/*        Y      - M-VECTOR TO WHICH AX WILL BE ADDED. */
/*        APNT   - INDEX VECTOR FOR A.  XA(I) POINTS TO THE */
/*                 FIRST NONZERO IN COLUMN I OF A. */
/*        Y      - ON OUTPUT, CONTAINS Y = Y + AX. */

/* *********************************************************************** */

template <typename Scalar, typename Index>
void smxpy1(
  Index m, Index n, Scalar* __restrict y, Index* __restrict apnt,
  Scalar* __restrict a) {
  /* System generated locals */
  Index i__3, i__4, i__5;
  Scalar z__2;

  /* Local variables */
  Index i, j, ii, jj;
  Scalar amult;

  /* Parameter adjustments */
  --y;
  --apnt;
  --a;

  /* Function Body */
  for (j = 1; j <= n; ++j) {
    jj = apnt[j];
    ii = apnt[j + 1] - m;
    z__2 = -a[jj];
    i__3 = ii;
    amult = z__2 * a[i__3];
    for (i = 1; i <= m; ++i) {
      i__3 = i;
      i__4 = i;
      i__5 = ii;
      z__2 = amult * a[i__5];
      y[i__3] = y[i__4] + z__2;
      ++ii;
    }
  }
} /* smxpy1 */

}}}  // namespace NgPeytonCpp::details::f2c

#endif  // NGPEYTONCPP_DETAILS_BLAS_H
