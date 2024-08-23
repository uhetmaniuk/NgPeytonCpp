/* igathr.f -- translated by f2c (version 20230428).
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
/* ******         IGATHR .... INTEGER GATHER OPERATION      ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS A STANDARD INTEGER GATHER */
/*               OPERATION. */

/*     INPUT PARAMETERS - */
/*        KLEN   - LENGTH OF THE LIST OF GLOBAL INDICES. */
/*        LINDX  - LIST OF GLOBAL INDICES. */
/*        INDMAP - INDEXED BY GLOBAL INDICES, IT CONTAINS THE */
/*                 REQUIRED RELATIVE INDICES. */

/*     OUTPUT PARAMETERS - */
/*        RELIND - LIST RELATIVE INDICES. */

/* *********************************************************************** */

/* Subroutine */ int igathr_(integer *klen, integer *lindx, integer *indmap, 
	integer *relind)
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

/* DIR$ IVDEP */
    /* Parameter adjustments */
    --relind;
    --indmap;
    --lindx;

    /* Function Body */
    i__1 = *klen;
    for (i__ = 1; i__ <= i__1; ++i__) {
	relind[i__] = indmap[lindx[i__]];
/* L100: */
    }
    return 0;
} /* igathr_ */

