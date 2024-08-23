/* mmpy2.f -- translated by f2c (version 20230428).
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
/* *************     MMPY2  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*       LOOP UNROLLING: LEVEL 2 */

/*   INPUT PARAMETERS - */
/*       M               -   NUMBER OF ROWS IN X AND IN Y. */
/*       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS */
/*                           IN A. */
/*       Q               -   NUMBER OF COLUMNS IN A AND Y. */
/*       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE */
/*                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO */
/*                           USED TO ACCESS THE ROWS OF A. */
/*       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A. */
/*       LDY             -   LENGTH OF FIRST COLUMN OF Y. */

/*   UPDATED PARAMETERS - */
/*       Y(*)            -   ON OUTPUT, Y = Y + AX. */

/* *********************************************************************** */

/* Subroutine */ int mmpy2_(integer *m, integer *n, integer *q, integer *xpnt,
	 doublereal *x, doublereal *y, integer *ldy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static doublereal a1, a2;
    static integer i1, i2, j1, j2, mm, iy, xcol, ycol, leny, remain, iylast, 
	    iystop, iystrt;


/* *********************************************************************** */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */


/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */



/* *********************************************************************** */

/*       ----------------------------------------------------------- */
/*       INITIAL OFFSETS, COLUMN LENGTHS, AND INDEX RANGE VARIABLES. */
/*       ----------------------------------------------------------- */
    /* Parameter adjustments */
    --y;
    --x;
    --xpnt;

    /* Function Body */
    remain = *n % 2 + 1;
    mm = *m;
    iylast = 0;
    leny = *ldy;

/*       ------------------------------------ */
/*       TO COMPUTE EACH COLUMN YCOL OF Y ... */
/*       ------------------------------------ */

    i__1 = *q;
    for (ycol = 1; ycol <= i__1; ++ycol) {

	iystrt = iylast + 1;
	iystop = iystrt + mm - 1;
	iylast += leny;

/*           -------------------------------------------------- */
/*           ... PERFORM THE APPROPRATE MATRIX VECTOR MULTIPLY: */
/*               X * A(*,YCOL) WITH LEVEL 2 LOOP-UNROLLING. */
/*           -------------------------------------------------- */
	switch (remain) {
	    case 1:  goto L200;
	    case 2:  goto L100;
	}

L100:
	j1 = xpnt[1];
	i1 = xpnt[2] - mm;
	a1 = -x[j1] * x[i1];
	i__2 = iystop;
	for (iy = iystrt; iy <= i__2; ++iy) {
	    y[iy] += a1 * x[i1];
	    ++i1;
/* L150: */
	}
	goto L200;

L200:
	i__2 = *n;
	for (xcol = remain; xcol <= i__2; xcol += 2) {
	    j1 = xpnt[xcol];
	    j2 = xpnt[xcol + 1];
	    i1 = j2 - mm;
	    i2 = xpnt[xcol + 2] - mm;
	    a1 = -x[j1] * x[i1];
	    a2 = -x[j2] * x[i2];
	    i__3 = iystop;
	    for (iy = iystrt; iy <= i__3; ++iy) {
		y[iy] = y[iy] + a1 * x[i1] + a2 * x[i2];
		++i1;
		++i2;
/* L300: */
	    }
/* L400: */
	}

	--mm;
	--leny;

/* L500: */
    }

    return 0;
} /* mmpy2_ */

