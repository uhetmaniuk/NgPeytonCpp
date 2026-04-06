#ifndef NGPEYTONCPP_DETAILS_SOLVE_H
#define NGPEYTONCPP_DETAILS_SOLVE_H

namespace NgPeytonCpp { namespace details { namespace f2c {

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********     BLKSLV ... BLOCK TRIANGULAR SOLUTIONS          ********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       GIVEN THE CHOLESKY FACTORIZATION OF A SPARSE SYMMETRIC */
/*       POSITIVE DEFINITE MATRIX, THIS SUBROUTINE PERFORMS THE */
/*       TRIANGULAR SOLUTION.  IT USES OUTPUT FROM BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   SUPERNODE PARTITION. */
/*       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE. */
/*       (XLNZ,LNZ)      -   CHOLESKY FACTOR. */

/*   UPDATED PARAMETERS: */
/*       RHS             -   ON INPUT, CONTAINS THE RIGHT HAND SIDE.  ON */
/*                           OUTPUT, CONTAINS THE SOLUTION. */

/* *********************************************************************** */

template <typename Scalar, typename Index>
void blkslv(
  Index nsuper, Index* __restrict xsuper, Index* __restrict xlindx,
  Index* __restrict lindx, Index* __restrict xlnz, Scalar* __restrict lnz,
  Scalar* __restrict rhs) {
  /* System generated locals */
  Index i__3, i__4, i__5, i__6;
  Scalar z__1, z__2;

  /* Local variables */
  Index i;
  Scalar t;
  Index ix, jcol, ipnt, jpnt, jsup, fjcol, ljcol, ixstop, ixstrt;

  /* Parameter adjustments */
  --rhs;
  --lnz;
  --xlnz;
  --lindx;
  --xlindx;
  --xsuper;

  /* Function Body */
  if (nsuper <= 0) {
    return;
  }

  /*       ------------------------ */
  /*       FORWARD SUBSTITUTION ... */
  /*       ------------------------ */
  fjcol = xsuper[1];
  for (jsup = 1; jsup <= nsuper; ++jsup) {
    ljcol = xsuper[jsup + 1] - 1;
    ixstrt = xlnz[fjcol];
    jpnt = xlindx[jsup];
    for (jcol = fjcol; jcol <= ljcol; ++jcol) {
      ixstop = xlnz[jcol + 1] - 1;
      t = rhs[jcol];
      ipnt = jpnt + 1;
      for (ix = ixstrt + 1; ix <= ixstop; ++ix) {
        i = lindx[ipnt];
        i__4 = i;
        i__5 = i;
        i__6 = ix;
        z__2 = t * lnz[i__6];
        z__1 = rhs[i__5] - z__2;
        rhs[i__4] = z__1;
        ++ipnt;
      }
      ixstrt = ixstop + 1;
      ++jpnt;
    }
    fjcol = ljcol + 1;
  }

  /*       ------------------ */
  /*       DIAGONAL SOLVE ... */
  /*       ------------------ */
  for (jcol = 1; jcol < xsuper[nsuper + 1]; ++jcol) {
    rhs[jcol] = rhs[jcol] / lnz[xlnz[jcol]];
  }

  /*       ------------------------- */
  /*       BACKWARD SUBSTITUTION ... */
  /*       ------------------------- */
  ljcol = xsuper[nsuper + 1] - 1;
  for (jsup = nsuper; jsup >= 1; --jsup) {
    fjcol = xsuper[jsup];
    ixstop = xlnz[ljcol + 1] - 1;
    jpnt = xlindx[jsup] + (ljcol - fjcol);
    for (jcol = ljcol; jcol >= fjcol; --jcol) {
      ixstrt = xlnz[jcol];
      ipnt = jpnt + 1;
      t = rhs[jcol];
      for (ix = ixstrt + 1; ix <= ixstop; ++ix) {
        i = lindx[ipnt];
        i__3 = ix;
        i__4 = i;
        z__2 = lnz[i__3] * rhs[i__4];
        z__1 = t - z__2;
        t = z__1;
        ++ipnt;
      }
      rhs[jcol] = t;
      ixstop = ixstrt - 1;
      --jpnt;
    }
    ljcol = fjcol - 1;
  }
} /* blkslv_ */

}}}  // namespace NgPeytonCpp::details::f2c

#endif  // NGPEYTONCPP_DETAILS_SOLVE_H
