/* etpost.f -- translated by f2c (version 20230428).
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
/* ***************     ETPOST ..... ETREE POSTORDERING     *************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (SEPT 17, 1986) */

/*   PURPOSE: */
/*       BASED ON THE BINARY REPRESENTATION (FIRST-SON,BROTHER) OF */
/*       THE ELIMINATION TREE, A POSTORDERING IS DETERMINED. THE */
/*       CORRESPONDING PARENT VECTOR IS ALSO MODIFIED TO REFLECT */
/*       THE REORDERING. */

/*   INPUT PARAMETERS: */
/*       ROOT            -   ROOT OF THE ELIMINATION TREE (USUALLY IT */
/*                           IS NEQNS). */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHR VECTOR. */

/*   UPDATED PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR. */

/*   OUTPUT PARAMETERS: */
/*       INVPOS          -   INVERSE PERMUTATION FOR THE POSTORDERING. */

/*   WORKING PARAMETERS: */
/*       STACK           -   THE STACK FOR POSTORDER TRAVERSAL OF THE */
/*                           TREE. */

/* *********************************************************************** */

/* Subroutine */ int etpost_(integer *root, integer *fson, integer *brothr, 
	integer *invpos, integer *parent, integer *stack)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer num, node, itop, ndpar, nunode;


/* *********************************************************************** */



/* *********************************************************************** */


/* *********************************************************************** */

    /* Parameter adjustments */
    --stack;
    --parent;
    --invpos;
    --brothr;
    --fson;

    /* Function Body */
    num = 0;
    itop = 0;
    node = *root;
/*       ------------------------------------------------------------- */
/*       TRAVERSE ALONG THE FIRST SONS POINTER AND PUSH THE TREE NODES */
/*       ALONG THE TRAVERSAL INTO THE STACK. */
/*       ------------------------------------------------------------- */
L100:
    ++itop;
    stack[itop] = node;
    node = fson[node];
    if (node > 0) {
	goto L100;
    }
/*           ---------------------------------------------------------- */
/*           IF POSSIBLE, POP A TREE NODE FROM THE STACK AND NUMBER IT. */
/*           ---------------------------------------------------------- */
L200:
    if (itop <= 0) {
	goto L300;
    }
    node = stack[itop];
    --itop;
    ++num;
    invpos[node] = num;
/*               ---------------------------------------------------- */
/*               THEN, TRAVERSE TO ITS YOUNGER BROTHER IF IT HAS ONE. */
/*               ---------------------------------------------------- */
    node = brothr[node];
    if (node <= 0) {
	goto L200;
    }
    goto L100;

L300:
/*       ------------------------------------------------------------ */
/*       DETERMINE THE NEW PARENT VECTOR OF THE POSTORDERING.  BROTHR */
/*       IS USED TEMPORARILY FOR THE NEW PARENT VECTOR. */
/*       ------------------------------------------------------------ */
    i__1 = num;
    for (node = 1; node <= i__1; ++node) {
	nunode = invpos[node];
	ndpar = parent[node];
	if (ndpar > 0) {
	    ndpar = invpos[ndpar];
	}
	brothr[nunode] = ndpar;
/* L400: */
    }

    i__1 = num;
    for (nunode = 1; nunode <= i__1; ++nunode) {
	parent[nunode] = brothr[nunode];
/* L500: */
    }

    return 0;
} /* etpost_ */

