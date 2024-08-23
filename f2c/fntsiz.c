/* fntsiz.f -- translated by f2c (version 20230428).
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
/* ******     FNTSIZ ..... COMPUTE WORK STORAGE SIZE FOR BLKFCT     ****** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES THE SIZE OF THE WORKING STORAGE */
/*       REQUIRED BY BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       SNODE           -   SUPERNODE MEMBERSHIP. */
/*       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE. */

/*   OUTPUT PARAMETERS: */
/*       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT. */

/* *********************************************************************** */

/* Subroutine */ int fntsiz_(integer *nsuper, integer *xsuper, integer *snode,
	 integer *xlindx, integer *lindx, integer *tmpsiz)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, iend, clen, ksup, bound, ncols, width, tsize, ibegin, 
	    length, cursup, nxtsup;


/* *********************************************************************** */



/* *********************************************************************** */

/*       RETURNS SIZE OF TEMP ARRAY USED BY BLKFCT FACTORIZATION ROUTINE. */
/*       NOTE THAT THE VALUE RETURNED IS AN ESTIMATE, THOUGH IT IS USUALLY */
/*       TIGHT. */

/*       ---------------------------------------- */
/*       COMPUTE SIZE OF TEMPORARY STORAGE VECTOR */
/*       NEEDED BY BLKFCT. */
/*       ---------------------------------------- */
    /* Parameter adjustments */
    --lindx;
    --xlindx;
    --snode;
    --xsuper;

    /* Function Body */
    *tmpsiz = 0;
    for (ksup = *nsuper; ksup >= 1; --ksup) {
	ncols = xsuper[ksup + 1] - xsuper[ksup];
	ibegin = xlindx[ksup] + ncols;
	iend = xlindx[ksup + 1] - 1;
	length = iend - ibegin + 1;
	bound = length * (length + 1) / 2;
	if (bound > *tmpsiz) {
	    cursup = snode[lindx[ibegin]];
	    clen = xlindx[cursup + 1] - xlindx[cursup];
	    width = 0;
	    i__1 = iend;
	    for (i__ = ibegin; i__ <= i__1; ++i__) {
		nxtsup = snode[lindx[i__]];
		if (nxtsup == cursup) {
		    ++width;
		    if (i__ == iend) {
			if (clen > length) {
			    tsize = length * width - (width - 1) * width / 2;
			    *tmpsiz = max(tsize,*tmpsiz);
			}
		    }
		} else {
		    if (clen > length) {
			tsize = length * width - (width - 1) * width / 2;
			*tmpsiz = max(tsize,*tmpsiz);
		    }
		    length -= width;
		    bound = length * (length + 1) / 2;
		    if (bound <= *tmpsiz) {
			goto L500;
		    }
		    width = 1;
		    cursup = nxtsup;
		    clen = xlindx[cursup + 1] - xlindx[cursup];
		}
/* L400: */
	    }
	}
L500:
	;
    }

    return 0;
} /* fntsiz_ */

