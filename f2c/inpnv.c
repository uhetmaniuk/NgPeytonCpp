/* inpnv.f -- translated by f2c (version 20230428).
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

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */

/*     ------------------------------------------------------ */
/*     INPUT NUMERICAL VALUES INTO SPARSE DATA STRUCTURES ... */
/*     ------------------------------------------------------ */

/* Subroutine */ int inpnv_(integer *neqns, integer *xadjf, 
	integer *adjf, doublereal *anzf, integer *perm, integer *invp, 
	integer *nsuper, integer *xsuper, integer *xlindx, integer *lindx, 
	integer *xlnz, doublereal *lnz, integer *iwsiz, integer *offset, 
	integer *iflag)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, ii, jlen, oldj, last, jsuper;

    /* Parameter adjustments */
    --offset;
    --lnz;
    --xlnz;
    --lindx;
    --xlindx;
    --xsuper;
    --invp;
    --perm;
    --anzf;
    --adjf;
    --xadjf;

    /* Function Body */
    *iflag = 0;
    if (*iwsiz < *neqns) {
	*iflag = -1;
	return 0;
    }

    i__1 = *nsuper;
    for (jsuper = 1; jsuper <= i__1; ++jsuper) {

/*           ---------------------------------------- */
/*           FOR EACH SUPERNODE, DO THE FOLLOWING ... */
/*           ---------------------------------------- */

/*           ----------------------------------------------- */
/*           FIRST GET OFFSET TO FACILITATE NUMERICAL INPUT. */
/*           ----------------------------------------------- */
	jlen = xlindx[jsuper + 1] - xlindx[jsuper];
	i__2 = xlindx[jsuper + 1] - 1;
	for (ii = xlindx[jsuper]; ii <= i__2; ++ii) {
	    i__ = lindx[ii];
	    --jlen;
	    offset[i__] = jlen;
/* L100: */
	}

	i__2 = xsuper[jsuper + 1] - 1;
	for (j = xsuper[jsuper]; j <= i__2; ++j) {
/*               ----------------------------------------- */
/*               FOR EACH COLUMN IN THE CURRENT SUPERNODE, */
/*               FIRST INITIALIZE THE DATA STRUCTURE. */
/*               ----------------------------------------- */
	    i__3 = xlnz[j + 1] - 1;
	    for (ii = xlnz[j]; ii <= i__3; ++ii) {
		lnz[ii] = 0.;
/* L200: */
	    }

/*               ----------------------------------- */
/*               NEXT INPUT THE INDIVIDUAL NONZEROS. */
/*               ----------------------------------- */
	    oldj = perm[j];
	    last = xlnz[j + 1] - 1;
	    i__3 = xadjf[oldj + 1] - 1;
	    for (ii = xadjf[oldj]; ii <= i__3; ++ii) {
		i__ = invp[adjf[ii]];
		if (i__ >= j) {
		    lnz[last - offset[i__]] = anzf[ii];
		}
/* L300: */
	    }
/* L400: */
	}

/* L500: */
    }
    return 0;
} /* inpnv_ */

