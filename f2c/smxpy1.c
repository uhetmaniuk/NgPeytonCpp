/* smxpy1.f -- translated by f2c (version 20230428).
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

/* Subroutine */ int smxpy1_(integer *m, integer *n, doublereal *y, integer *
	apnt, doublereal *a)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, ii, jj;
    static doublereal amult;


/* *********************************************************************** */

/*     ----------- */
/*     PARAMETERS. */
/*     ----------- */




/*     ---------------- */
/*     LOCAL VARIABLES. */
/*     ---------------- */



/* *********************************************************************** */

    /* Parameter adjustments */
    --y;
    --apnt;
    --a;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jj = apnt[j];
	ii = apnt[j + 1] - *m;
	amult = -a[jj] * a[ii];
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__] += amult * a[ii];
	    ++ii;
/* L100: */
	}
/* L200: */
    }
    return 0;
} /* smxpy1_ */

