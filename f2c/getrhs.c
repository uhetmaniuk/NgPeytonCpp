/* getrhs.f -- translated by f2c (version 20230428).
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

/*     ------------------------------------ */
/*     CONSTRUCT RIGHT HAND SIDE VECTOR ... */
/*     ------------------------------------ */

/* Subroutine */ int getrhs_(integer *n, integer *xadjf, integer *adjf, 
	doublereal *anzf, doublereal *sol, doublereal *rhs)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j;
    static doublereal t;
    static integer ii;




/*       --------------- */
/*       INITIALIZATION. */
/*       --------------- */
    /* Parameter adjustments */
    --rhs;
    --sol;
    --anzf;
    --adjf;
    --xadjf;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	rhs[j] = 0.f;
/* L100: */
    }

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*           ------------------- */
/*           FOR EACH COLUMN ... */
/*           ------------------- */
	t = sol[j];
	i__2 = xadjf[j + 1] - 1;
	for (ii = xadjf[j]; ii <= i__2; ++ii) {
	    rhs[adjf[ii]] += t * anzf[ii];
/* L200: */
	}
/* L300: */
    }

    return 0;
} /* getrhs_ */

