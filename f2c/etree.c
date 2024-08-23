/* etree.f -- translated by f2c (version 20230428).
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
/* ****************     ETREE ..... ELIMINATION TREE     ***************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO DETERMINE THE ELIMINATION TREE FROM A GIVEN ORDERING AND */
/*       THE ADJACENCY STRUCTURE.  THE PARENT VECTOR IS RETURNED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE. */
/*       (PERM,INVP)     -   PERMUTATION AND INVERSE PERMUTATION VECTORS */

/*   OUTPUT PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE. */

/*   WORKING PARAMETERS: */
/*       ANCSTR          -   THE ANCESTOR VECTOR. */

/* *********************************************************************** */

/* Subroutine */ int etree_(integer *neqns, integer *xadj, integer *adjncy, 
	integer *perm, integer *invp, integer *parent, integer *ancstr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, nbr, node, next, jstop, jstrt;


/* *********************************************************************** */



/* *********************************************************************** */


/* *********************************************************************** */

    /* Parameter adjustments */
    --ancstr;
    --parent;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    if (*neqns <= 0) {
	return 0;
    }

    i__1 = *neqns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	parent[i__] = 0;
	ancstr[i__] = 0;
	node = perm[i__];

	jstrt = xadj[node];
	jstop = xadj[node + 1] - 1;
	if (jstrt <= jstop) {
	    i__2 = jstop;
	    for (j = jstrt; j <= i__2; ++j) {
		nbr = adjncy[j];
		nbr = invp[nbr];
		if (nbr < i__) {
/*                       ------------------------------------------- */
/*                       FOR EACH NBR, FIND THE ROOT OF ITS CURRENT */
/*                       ELIMINATION TREE.  PERFORM PATH COMPRESSION */
/*                       AS THE SUBTREE IS TRAVERSED. */
/*                       ------------------------------------------- */
L100:
		    if (ancstr[nbr] == i__) {
			goto L300;
		    }
		    if (ancstr[nbr] > 0) {
			next = ancstr[nbr];
			ancstr[nbr] = i__;
			nbr = next;
			goto L100;
		    }
/*                       -------------------------------------------- */
/*                       NOW, NBR IS THE ROOT OF THE SUBTREE.  MAKE I */
/*                       THE PARENT NODE OF THIS ROOT. */
/*                       -------------------------------------------- */
		    parent[nbr] = i__;
		    ancstr[nbr] = i__;
		}
L300:
		;
	    }
	}
/* L400: */
    }

    return 0;
} /* etree_ */

