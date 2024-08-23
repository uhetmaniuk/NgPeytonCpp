/* invinv.f -- translated by f2c (version 20230428).
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
/*   Authors:        Joseph W.H. Liu */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ***********     INVINV ..... CONCATENATION OF TWO INVP     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO PERFORM THE MAPPING OF */
/*           ORIGINAL-INVP --> INTERMEDIATE-INVP --> NEW INVP */
/*       AND THE RESULTING ORDERING REPLACES INVP.  THE NEW PERMUTATION */
/*       VECTOR PERM IS ALSO COMPUTED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       INVP2           -   THE SECOND INVERSE PERMUTATION VECTOR. */

/*   UPDATED PARAMETERS: */
/*       INVP            -   THE FIRST INVERSE PERMUTATION VECTOR.  ON */
/*                           OUTPUT, IT CONTAINS THE NEW INVERSE */
/*                           PERMUTATION. */

/*   OUTPUT PARAMETER: */
/*       PERM            -   NEW PERMUTATION VECTOR (CAN BE THE SAME AS */
/*                           INVP2). */

/* *********************************************************************** */

/* Subroutine */ int invinv_(integer *neqns, integer *invp, integer *invp2, 
	integer *perm)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, node, interm;


/* *********************************************************************** */



/* *********************************************************************** */


/* *********************************************************************** */

    /* Parameter adjustments */
    --perm;
    --invp2;
    --invp;

    /* Function Body */
    i__1 = *neqns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	interm = invp[i__];
	invp[i__] = invp2[interm];
/* L100: */
    }

    i__1 = *neqns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	node = invp[i__];
	perm[node] = i__;
/* L200: */
    }

    return 0;
} /* invinv_ */

