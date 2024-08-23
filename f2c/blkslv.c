/* blkslv.f -- translated by f2c (version 20230428).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

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

/* Subroutine */ int blkslv_(integer *nsuper, integer *xsuper, integer *
	xlindx, integer *lindx, integer *xlnz, doublereal *lnz, doublereal *
	rhs)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__;
    static doublereal t;
    static integer ix, jcol, ipnt, jpnt, jsup, fjcol, ljcol, ixstop, ixstrt;


/* *********************************************************************** */


/* *********************************************************************** */


/* *********************************************************************** */

    /* Parameter adjustments */
    --rhs;
    --lnz;
    --xlnz;
    --lindx;
    --xlindx;
    --xsuper;

    /* Function Body */
    if (*nsuper <= 0) {
	return 0;
    }

/*       ------------------------ */
/*       FORWARD SUBSTITUTION ... */
/*       ------------------------ */
    fjcol = xsuper[1];
    i__1 = *nsuper;
    for (jsup = 1; jsup <= i__1; ++jsup) {
	ljcol = xsuper[jsup + 1] - 1;
	ixstrt = xlnz[fjcol];
	jpnt = xlindx[jsup];
	i__2 = ljcol;
	for (jcol = fjcol; jcol <= i__2; ++jcol) {
	    ixstop = xlnz[jcol + 1] - 1;
	    t = rhs[jcol];
	    ipnt = jpnt + 1;
/* DIR$           IVDEP */
	    i__3 = ixstop;
	    for (ix = ixstrt + 1; ix <= i__3; ++ix) {
		i__ = lindx[ipnt];
		rhs[i__] -= t * lnz[ix];
		++ipnt;
/* L100: */
	    }
	    ixstrt = ixstop + 1;
	    ++jpnt;
/* L200: */
	}
	fjcol = ljcol + 1;
/* L300: */
    }

/*       ------------------ */
/*       DIAGONAL SOLVE ... */
/*       ------------------ */
    i__1 = xsuper[*nsuper + 1] - 1;
    for (jcol = 1; jcol <= i__1; ++jcol) {
	rhs[jcol] /= lnz[xlnz[jcol]];
/* L400: */
    }

/*       ------------------------- */
/*       BACKWARD SUBSTITUTION ... */
/*       ------------------------- */
    ljcol = xsuper[*nsuper + 1] - 1;
    for (jsup = *nsuper; jsup >= 1; --jsup) {
	fjcol = xsuper[jsup];
	ixstop = xlnz[ljcol + 1] - 1;
	jpnt = xlindx[jsup] + (ljcol - fjcol);
	i__1 = fjcol;
	for (jcol = ljcol; jcol >= i__1; --jcol) {
	    ixstrt = xlnz[jcol];
	    ipnt = jpnt + 1;
	    t = rhs[jcol];
/* DIR$           IVDEP */
	    i__2 = ixstop;
	    for (ix = ixstrt + 1; ix <= i__2; ++ix) {
		i__ = lindx[ipnt];
		t -= lnz[ix] * rhs[i__];
		++ipnt;
/* L500: */
	    }
	    rhs[jcol] = t;
	    ixstop = ixstrt - 1;
	    --jpnt;
/* L600: */
	}
	ljcol = fjcol - 1;
/* L700: */
    }

    return 0;
} /* blkslv_ */

