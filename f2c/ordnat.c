/* ordnat.f -- translated by f2c (version 20230428).
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
/* ****     ORDNAT ..... NATURAL ORDERING                     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE RECORDS THE INITIAL ORDERING IN THE */
/*               ORDERING VECTORS PERM AND INVP. */

/*     INPUT PARAMETERS - */
/*        NEQNS  - NUMBER OF EQUATIONS. */

/*     OUTPUT PARAMETERS - */
/*        PERM   - THE "NATURAL" ORDERING; I.E., THE INITIAL */
/*                 ORDERING. */
/*        INVP   - THE INVERSE OF PERM. */
/*        SFIFLG - SFIFLG=.F. MEANS SKIP SYMBOLIC FACTORIZATION */
/*                 INITIALIZATION (SFINIT), SFIFLG=.T. MEANS EXECUTE */
/*                 SFINIT. */

/* *********************************************************************** */

/* Subroutine */ int ordnat_(integer *neqns, integer *perm, integer *invp, 
	logical *sfiflg)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/* *********************************************************************** */



/* *********************************************************************** */

    /* Parameter adjustments */
    --invp;
    --perm;

    /* Function Body */
    i__1 = *neqns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	perm[i__] = i__;
	invp[i__] = i__;
/* L700: */
    }
    *sfiflg = TRUE_;
    return 0;

} /* ordnat_ */

