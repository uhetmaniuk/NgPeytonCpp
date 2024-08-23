/* fsup2.f -- translated by f2c (version 20230428).
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
/* ****************    FSUP2  ..... FIND SUPERNODES #2   ***************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE IS THE SECOND OF TWO ROUTINES FOR FINDING A */
/*       MAXIMAL SUPERNODE PARTITION.  IT'S SOLE PURPOSE IS TO */
/*       CONSTRUCT THE NEEDED VECTOR OF LENGTH NSUPER: XSUPER(*).  THE */
/*       FIRST ROUTINE FSUP1 COMPUTES THE NUMBER OF SUPERNODES AND THE */
/*       SUPERNODE MEMBERSHIP VECTOR SNODE(*), WHICH IS OF LENGTH NEQNS. */


/*   ASSUMPTIONS: */
/*       THIS ROUTINE ASSUMES A POSTORDERING OF THE ELIMINATION TREE.  IT */
/*       ALSO ASSUMES THAT THE OUTPUT FROM FSUP1 IS AVAILABLE. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS). */
/*       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           ELIMINATION TREE OF THE POSTORDERED MATRIX. */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */

/*   OUTPUT PARAMETERS: */
/*       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE */
/*                           SUPERNODE PARTITIONING. */

/*   FIRST CREATED ON    JANUARY 18, 1992. */
/*   LAST UPDATED ON     NOVEMEBER 22, 1994. */

/* *********************************************************************** */

/* Subroutine */ int fsup2_(integer *neqns, integer *nsuper, integer *etpar, 
	integer *snode, integer *xsuper)
{
    static integer kcol, ksup, lstsup;


/* *********************************************************************** */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* *********************************************************************** */

/*       ------------------------------------------------- */
/*       COMPUTE THE SUPERNODE PARTITION VECTOR XSUPER(*). */
/*       ------------------------------------------------- */
    /* Parameter adjustments */
    --xsuper;
    --snode;
    --etpar;

    /* Function Body */
    lstsup = *nsuper + 1;
    for (kcol = *neqns; kcol >= 1; --kcol) {
	ksup = snode[kcol];
	if (ksup != lstsup) {
	    xsuper[lstsup] = kcol + 1;
	}
	lstsup = ksup;
/* L100: */
    }
    xsuper[1] = 1;

    return 0;
} /* fsup2_ */

