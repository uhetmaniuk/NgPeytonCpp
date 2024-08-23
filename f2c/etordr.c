/* etordr.f -- translated by f2c (version 20230428).
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
/* **********     ETORDR ..... ELIMINATION TREE REORDERING     *********** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   WRITTEN BY JOSEPH LIU (JUL 17, 1985) */

/*   PURPOSE: */
/*       TO DETERMINE AN EQUIVALENT REORDERING BASED ON THE STRUCTURE OF */
/*       THE ELIMINATION TREE.  A POSTORDERING OF THE GIVEN ELIMINATION */
/*       TREE IS RETURNED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE. */

/*   UPDATED PARAMETERS: */
/*       (PERM,INVP)     -   ON INPUT, THE GIVEN PERM AND INVERSE PERM */
/*                           VECTORS.  ON OUTPUT, THE NEW PERM AND */
/*                           INVERSE PERM VECTORS OF THE EQUIVALENT */
/*                           ORDERING. */

/*   OUTPUT PARAMETERS: */
/*       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE */
/*                           ASSOCIATED WITH THE NEW ORDERING. */

/*   WORKING PARAMETERS: */
/*       FSON            -   THE FIRST SON VECTOR. */
/*       BROTHR          -   THE BROTHER VECTOR. */
/*       INVPOS          -   THE INVERSE PERM VECTOR FOR THE */
/*                           POSTORDERING. */

/*   PROGRAM SUBROUTINES: */
/*       BETREE, ETPOST, ETREE , INVINV. */

/* *********************************************************************** */

/* Subroutine */ int etordr_(integer *neqns, integer *xadj, integer *adjncy, 
	integer *perm, integer *invp, integer *parent, integer *fson, integer 
	*brothr, integer *invpos)
{
    extern /* Subroutine */ int etree_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), betree_(integer *, 
	    integer *, integer *, integer *), invinv_(integer *, integer *, 
	    integer *, integer *), etpost_(integer *, integer *, integer *, 
	    integer *, integer *, integer *);


/* *********************************************************************** */



/* *********************************************************************** */

/*       ----------------------------- */
/*       COMPUTE THE ELIMINATION TREE. */
/*       ----------------------------- */
    /* Parameter adjustments */
    --invpos;
    --brothr;
    --fson;
    --parent;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    etree_(neqns, &xadj[1], &adjncy[1], &perm[1], &invp[1], &parent[1], &
	    invpos[1]);

/*       -------------------------------------------------------- */
/*       COMPUTE A BINARY REPRESENTATION OF THE ELIMINATION TREE. */
/*       -------------------------------------------------------- */
    betree_(neqns, &parent[1], &fson[1], &brothr[1]);

/*       ------------------------------- */
/*       POSTORDER THE ELIMINATION TREE. */
/*       ------------------------------- */
    etpost_(neqns, &fson[1], &brothr[1], &invpos[1], &parent[1], &perm[1]);

/*       -------------------------------------------------------- */
/*       COMPOSE THE ORIGINAL ORDERING WITH THE NEW POSTORDERING. */
/*       -------------------------------------------------------- */
    invinv_(neqns, &invp[1], &invpos[1], &perm[1]);

    return 0;
} /* etordr_ */

