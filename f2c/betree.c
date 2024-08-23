/* betree.f -- translated by f2c (version 20230428).
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
/* ******     BETREE ..... BINARY TREE REPRESENTATION OF ETREE     ******* */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO DETERMINE THE BINARY TREE REPRESENTATION OF THE ELIMINATION */
/*       TREE GIVEN BY THE PARENT VECTOR.  THE RETURNED REPRESENTATION */
/*       WILL BE GIVEN BY THE FIRST-SON AND BROTHER VECTORS.  THE ROOT */
/*       OF THE BINARY TREE IS ALWAYS NEQNS. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE. */
/*                           IT IS ASSUMED THAT PARENT(I) > I EXCEPT OF */
/*                           THE ROOTS. */

/*   OUTPUT PARAMETERS: */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHER VECTOR. */

/* *********************************************************************** */

/* Subroutine */ int betree_(integer *neqns, integer *parent, integer *fson, 
	integer *brothr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer node, ndpar, lroot;


/* *********************************************************************** */



/* *********************************************************************** */


/* *********************************************************************** */

    /* Parameter adjustments */
    --brothr;
    --fson;
    --parent;

    /* Function Body */
    if (*neqns <= 0) {
	return 0;
    }

    i__1 = *neqns;
    for (node = 1; node <= i__1; ++node) {
	fson[node] = 0;
	brothr[node] = 0;
/* L100: */
    }
    lroot = *neqns;
/*       ------------------------------------------------------------ */
/*       FOR EACH NODE := NEQNS-1 STEP -1 DOWNTO 1, DO THE FOLLOWING. */
/*       ------------------------------------------------------------ */
    if (*neqns <= 1) {
	return 0;
    }
    for (node = *neqns - 1; node >= 1; --node) {
	ndpar = parent[node];
	if (ndpar <= 0 || ndpar == node) {
/*               ------------------------------------------------- */
/*               NODE HAS NO PARENT.  GIVEN STRUCTURE IS A FOREST. */
/*               SET NODE TO BE ONE OF THE ROOTS OF THE TREES. */
/*               ------------------------------------------------- */
	    brothr[lroot] = node;
	    lroot = node;
	} else {
/*               ------------------------------------------- */
/*               OTHERWISE, BECOMES FIRST SON OF ITS PARENT. */
/*               ------------------------------------------- */
	    brothr[node] = fson[ndpar];
	    fson[ndpar] = node;
	}
/* L300: */
    }
    brothr[lroot] = 0;

    return 0;
} /* betree_ */

