/* dscal.f -- translated by f2c (version 20230428).
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

/* Subroutine */ int dscal_(integer *n, doublereal *a, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/* *********************************************************************** */

/*     ----------- */
/*     PARAMETERS. */
/*     ----------- */

/*     ---------------- */
/*     LOCAL VARIABLES. */
/*     ---------------- */

/* *********************************************************************** */

    /* Parameter adjustments */
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = *a * x[i__];
/* L100: */
    }
    return 0;
} /* dscal_ */

