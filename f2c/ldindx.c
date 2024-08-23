/* ldindx.f -- translated by f2c (version 20230428).
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
/* ******         LDINDX .... LOAD INDEX VECTOR             ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE COMPUTES THE SECOND INDEX VECTOR */
/*               USED TO IMPLEMENT THE DOUBLY-INDIRECT SAXPY-LIKE */
/*               LOOPS THAT ALLOW US TO ACCUMULATE UPDATE */
/*               COLUMNS DIRECTLY INTO FACTOR STORAGE. */

/*     INPUT PARAMETERS - */
/*        JLEN   - LENGTH OF THE FIRST COLUMN OF THE SUPERNODE, */
/*                 INCLUDING THE DIAGONAL ENTRY. */
/*        LINDX  - THE OFF-DIAGONAL ROW INDICES OF THE SUPERNODE, */
/*                 I.E., THE ROW INDICES OF THE NONZERO ENTRIES */
/*                 LYING BELOW THE DIAGONAL ENTRY OF THE FIRST */
/*                 COLUMN OF THE SUPERNODE. */

/*     OUTPUT PARAMETERS - */
/*        INDMAP - THIS INDEX VECTOR MAPS EVERY GLOBAL ROW INDEX */
/*                 OF NONZERO ENTRIES IN THE FIRST COLUMN OF THE */
/*                 SUPERNODE TO ITS POSITION IN THE INDEX LIST */
/*                 RELATIVE TO THE LAST INDEX IN THE LIST.  MORE */
/*                 PRECISELY, IT GIVES THE DISTANCE OF EACH INDEX */
/*                 FROM THE LAST INDEX IN THE LIST. */

/* *********************************************************************** */

/* Subroutine */ int ldindx_(integer *jlen, integer *lindx, integer *indmap)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, jsub, curlen;


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
    --indmap;
    --lindx;

    /* Function Body */
    curlen = *jlen;
    i__1 = *jlen;
    for (j = 1; j <= i__1; ++j) {
	jsub = lindx[j];
	--curlen;
	indmap[jsub] = curlen;
/* L200: */
    }
    return 0;
} /* ldindx_ */

