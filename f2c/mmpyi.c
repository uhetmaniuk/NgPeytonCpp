/* mmpyi.f -- translated by f2c (version 20230428).
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
/* *************     MMPYI  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*       MATRIX X HAS ONLY 1 COLUMN. */

/*   INPUT PARAMETERS - */
/*       M               -   NUMBER OF ROWS IN X AND IN Y. */
/*       Q               -   NUMBER OF COLUMNS IN A AND Y. */
/*       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE */
/*                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO */
/*                           USED TO ACCESS THE ROWS OF A. */
/*       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A. */
/*       D               -   DIAGONAL ENTRY IN COLUMN X BACK IN CHOLESKY. */
/*       IY(*)           -   IY(COL) POINTS TO THE BEGINNING OF COLUMN */
/*       RELIND(*)       -   RELATIVE INDICES. */

/*   UPDATED PARAMETERS - */
/*       Y(*)            -   ON OUTPUT, Y = Y + AX. */

/* *********************************************************************** */

/* Subroutine */ int mmpyi_(integer *m, integer *q, integer *xpnt, doublereal 
	*x, doublereal *d__, integer *iy, doublereal *y, integer *relind)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal a;
    static integer i__, k, col, isub, ylast;


/* *********************************************************************** */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */


/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */


/* *********************************************************************** */

    /* Parameter adjustments */
    --relind;
    --y;
    --iy;
    --x;
    --xpnt;

    /* Function Body */
    i__1 = *q;
    for (k = 1; k <= i__1; ++k) {
	col = xpnt[k];
	ylast = iy[col + 1] - 1;
	a = -(*d__) * x[k];
/* DIR$   IVDEP */
	i__2 = *m;
	for (i__ = k; i__ <= i__2; ++i__) {
	    isub = xpnt[i__];
	    isub = ylast - relind[isub];
	    y[isub] += a * x[i__];
/* L100: */
	}
/* L200: */
    }
    return 0;

} /* mmpyi_ */

