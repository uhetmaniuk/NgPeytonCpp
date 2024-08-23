/* smxpy4.f -- translated by f2c (version 20230428).
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
/* ******     SMXPY4 .... MATRIX-VECTOR MULTIPLY            ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS A MATRIX-VECTOR MULTIPLY, */
/*               Y = Y + AX, ASSUMING DATA STRUCTURES USED IN */
/*               RECENTLY DEVELOPED SPARSE CHOLESKY CODES.  THE */
/*               '4' SIGNIFIES LEVEL 4 LOOP UNROLLING. */

/*     INPUT PARAMETERS - */
/*        M      - NUMBER OF ROWS. */
/*        N      - NUMBER OF COLUMNS. */
/*        Y      - M-VECTOR TO WHICH AX WILL BE ADDED. */
/*        APNT   - INDEX VECTOR FOR A.  XA(I) POINTS TO THE */
/*                 FIRST NONZERO IN COLUMN I OF A. */
/*        Y      - ON OUTPUT, CONTAINS Y = Y + AX. */

/* *********************************************************************** */

/* Subroutine */ int smxpy4_(integer *m, integer *n, doublereal *y, integer *
	apnt, doublereal *a)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal a1, a2, a3, a4;
    static integer i1, i2, i3, i4, j1, j2, j3, j4, remain;


/* *********************************************************************** */

/*     ----------- */
/*     PARAMETERS. */
/*     ----------- */





/*     ---------------- */
/*     LOCAL VARIABLES. */
/*     ---------------- */



/* *********************************************************************** */

    /* Parameter adjustments */
    --a;
    --apnt;
    --y;

    /* Function Body */
    remain = *n % 4;

    switch (remain + 1) {
	case 1:  goto L2000;
	case 2:  goto L100;
	case 3:  goto L200;
	case 4:  goto L300;
    }

L100:
    j1 = apnt[1];
    i1 = apnt[2] - *m;
    a1 = -a[j1] * a[i1];
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] += a1 * a[i1];
	++i1;
/* L150: */
    }
    goto L2000;

L200:
    j1 = apnt[1];
    j2 = apnt[2];
    i1 = j2 - *m;
    i2 = apnt[3] - *m;
    a1 = -a[j1] * a[i1];
    a2 = -a[j2] * a[i2];
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = y[i__] + a1 * a[i1] + a2 * a[i2];
	++i1;
	++i2;
/* L250: */
    }
    goto L2000;

L300:
    j1 = apnt[1];
    j2 = apnt[2];
    j3 = apnt[3];
    i1 = j2 - *m;
    i2 = j3 - *m;
    i3 = apnt[4] - *m;
    a1 = -a[j1] * a[i1];
    a2 = -a[j2] * a[i2];
    a3 = -a[j3] * a[i3];
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = y[i__] + a1 * a[i1] + a2 * a[i2] + a3 * a[i3];
	++i1;
	++i2;
	++i3;
/* L350: */
    }
    goto L2000;

L2000:
    i__1 = *n;
    for (j = remain + 1; j <= i__1; j += 4) {
	j1 = apnt[j];
	j2 = apnt[j + 1];
	j3 = apnt[j + 2];
	j4 = apnt[j + 3];
	i1 = j2 - *m;
	i2 = j3 - *m;
	i3 = j4 - *m;
	i4 = apnt[j + 4] - *m;
	a1 = -a[j1] * a[i1];
	a2 = -a[j2] * a[i2];
	a3 = -a[j3] * a[i3];
	a4 = -a[j4] * a[i4];
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__] = y[i__] + a1 * a[i1] + a2 * a[i2] + a3 * a[i3] + a4 * a[
		    i4];
	    ++i1;
	    ++i2;
	    ++i3;
	    ++i4;
/* L3000: */
	}
/* L4000: */
    }

    return 0;
} /* smxpy4_ */

