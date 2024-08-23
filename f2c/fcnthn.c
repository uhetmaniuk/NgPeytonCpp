/* fcnthn.f -- translated by f2c (version 20230428).
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
/* **************     FCNTHN  ..... FIND NONZERO COUNTS    *************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES THE ROW COUNTS AND COLUMN COUNTS IN */
/*       THE CHOLESKY FACTOR.  IT USES A DISJOINT SET UNION ALGORITHM. */

/*       TECHNIQUES: */
/*       1) SUPERNODE DETECTION. */
/*       2) PATH HALVING. */
/*       3) NO UNION BY RANK. */

/*   NOTES: */
/*       1) ASSUMES A POSTORDERING OF THE ELIMINATION TREE. */
/*       2) IT ASSUMES NO DIAGONAL ENTRIES IN THE ADJACENCY STRUCTURE, */
/*          I.E., NO SELF-LOOPS. */

/*   INPUT PARAMETERS: */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS. */
/*       (I) ADJLEN      -   LENGTH OF ADJACENCY STRUCTURE. */
/*       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           TO THE ADJACENCY STRUCTURE. */
/*       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING */
/*                           THE ADJACENCY STRUCTURE, EXCLUDING THE */
/*                           DIAGONAL ENTRIES. */
/*       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           POSTORDERING. */
/*       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           INVERSE OF THE POSTORDERING. */
/*       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                           ELIMINATION TREE OF THE POSTORDERED MATRIX. */

/*   OUTPUT PARAMETERS: */
/*       (I) ROWCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH ROW OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) NLNZ        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING */
/*                           THE DIAGONAL ENTRIES. */

/*   WORK PARAMETERS: */
/*       (I) SET(*)      -   ARRAY OF LENGTH NEQNS USED TO MAINTAIN THE */
/*                           DISJOINT SETS (I.E., SUBTREES). */
/*       (I) PRVLF(*)    -   ARRAY OF LENGTH NEQNS USED TO RECORD THE */
/*                           PREVIOUS LEAF OF EACH ROW SUBTREE. */
/*       (I) LEVEL(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE LEVEL */
/*                           (DISTANCE FROM THE ROOT). */
/*       (I) WEIGHT(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING WEIGHTS */
/*                           USED TO COMPUTE COLUMN COUNTS. */
/*       (I) FDESC(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE */
/*                           FIRST (I.E., LOWEST-NUMBERED) DESCENDANT. */
/*       (I) NCHILD(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE */
/*                           NUMBER OF CHILDREN. */
/*       (I) PRVNBR(*)   -   ARRAY OF LENGTH NEQNS USED TO RECORD THE */
/*                           PREVIOUS ``LOWER NEIGHBOR'' OF EACH NODE. */

/*   FIRST CREATED ON    APRIL 12, 1990. */
/*   LAST UPDATED ON     JANUARY 12, 1995. */
/*   LAST UPDATED ON     OCTOBER 20, 1997. */

/* *********************************************************************** */

/* Subroutine */ int fcnthn_(integer *neqns, integer *adjlen, integer *xadj, 
	integer *adjncy, integer *perm, integer *invp, integer *etpar, 
	integer *rowcnt, integer *colcnt, integer *nlnz, integer *set, 
	integer *prvlf, integer *level, integer *weight, integer *fdesc, 
	integer *nchild, integer *prvnbr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, lca, temp, xsup, last1, last2, lflag, pleaf, hinbr, 
	    jstop, jstrt, ifdesc, oldnbr, parent, lownbr;


/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* *********************************************************************** */

/*       -------------------------------------------------- */
/*       COMPUTE LEVEL(*), FDESC(*), NCHILD(*). */
/*       INITIALIZE XSUP, ROWCNT(*), COLCNT(*), */
/*                  SET(*), PRVLF(*), WEIGHT(*), PRVNBR(*). */
/*       -------------------------------------------------- */
    /* Parameter adjustments */
    --prvnbr;
    --prvlf;
    --set;
    --colcnt;
    --rowcnt;
    --etpar;
    --invp;
    --perm;
    --adjncy;
    --xadj;

    /* Function Body */
    level[0] = 0;
    for (k = *neqns; k >= 1; --k) {
	rowcnt[k] = 1;
	colcnt[k] = 0;
	set[k] = k;
	prvlf[k] = 0;
	level[k] = level[etpar[k]] + 1;
	weight[k] = 1;
	fdesc[k] = k;
	nchild[k] = 0;
	prvnbr[k] = 0;
/* L100: */
    }
    nchild[0] = 0;
    fdesc[0] = 0;
    i__1 = *neqns;
    for (k = 1; k <= i__1; ++k) {
	parent = etpar[k];
	weight[parent] = 0;
	++nchild[parent];
	ifdesc = fdesc[k];
	if (ifdesc < fdesc[parent]) {
	    fdesc[parent] = ifdesc;
	}
/* L200: */
    }
/*       ------------------------------------ */
/*       FOR EACH ``LOW NEIGHBOR'' LOWNBR ... */
/*       ------------------------------------ */
    i__1 = *neqns;
    for (lownbr = 1; lownbr <= i__1; ++lownbr) {
	lflag = 0;
	ifdesc = fdesc[lownbr];
	oldnbr = perm[lownbr];
	jstrt = xadj[oldnbr];
	jstop = xadj[oldnbr + 1] - 1;
/*           -------------------------------------------- */
/*           ISOLATED VERTEX IS LEAF IN OWN LEAF SUBTREE. */
/*           -------------------------------------------- */
	if (jstrt > jstop) {
	    lflag = 1;
	}
/*           ----------------------------------------------- */
/*           FOR EACH ``HIGH NEIGHBOR'', HINBR OF LOWNBR ... */
/*           ----------------------------------------------- */
	i__2 = jstop;
	for (j = jstrt; j <= i__2; ++j) {
	    hinbr = invp[adjncy[j]];
	    if (hinbr > lownbr) {
		if (ifdesc > prvnbr[hinbr]) {
/*                       ------------------------- */
/*                       INCREMENT WEIGHT(LOWNBR). */
/*                       ------------------------- */
		    ++weight[lownbr];
		    pleaf = prvlf[hinbr];
/*                       ----------------------------------------- */
/*                       IF HINBR HAS NO PREVIOUS ``LOW NEIGHBOR'' */
/*                       THEN ... */
/*                       ----------------------------------------- */
		    if (pleaf == 0) {
/*                           ----------------------------------------- */
/*                           ... ACCUMULATE LOWNBR-->HINBR PATH LENGTH */
/*                               IN ROWCNT(HINBR). */
/*                           ----------------------------------------- */
			rowcnt[hinbr] = rowcnt[hinbr] + level[lownbr] - level[
				hinbr];
		    } else {
/*                           ----------------------------------------- */
/*                           ... OTHERWISE, LCA <-- FIND(PLEAF), WHICH */
/*                               IS THE LEAST COMMON ANCESTOR OF PLEAF */
/*                               AND LOWNBR. */
/*                               (PATH HALVING.) */
/*                           ----------------------------------------- */
			last1 = pleaf;
			last2 = set[last1];
			lca = set[last2];
L300:
			if (lca != last2) {
			    set[last1] = lca;
			    last1 = lca;
			    last2 = set[last1];
			    lca = set[last2];
			    goto L300;
			}
/*                           ------------------------------------- */
/*                           ACCUMULATE PLEAF-->LCA PATH LENGTH IN */
/*                           ROWCNT(HINBR). */
/*                           DECREMENT WEIGHT(LCA). */
/*                           ------------------------------------- */
			rowcnt[hinbr] = rowcnt[hinbr] + level[lownbr] - level[
				lca];
			--weight[lca];
		    }
/*                       ---------------------------------------------- */
/*                       LOWNBR NOW BECOMES ``PREVIOUS LEAF'' OF HINBR. */
/*                       ---------------------------------------------- */
		    prvlf[hinbr] = lownbr;
		    lflag = 1;
		}
/*                   -------------------------------------------------- */
/*                   LOWNBR NOW BECOMES ``PREVIOUS NEIGHBOR'' OF HINBR. */
/*                   -------------------------------------------------- */
		prvnbr[hinbr] = lownbr;
	    }
/* L500: */
	}
/*           ---------------------------------------------------- */
/*           DECREMENT WEIGHT ( PARENT(LOWNBR) ). */
/*           SET ( P(LOWNBR) ) <-- SET ( P(LOWNBR) ) + SET(XSUP). */
/*           ---------------------------------------------------- */
	parent = etpar[lownbr];
	--weight[parent];
	if (lflag == 1 || nchild[lownbr] >= 2) {
	    xsup = lownbr;
	}
	set[xsup] = parent;
/* L600: */
    }
/*       --------------------------------------------------------- */
/*       USE WEIGHTS TO COMPUTE COLUMN (AND TOTAL) NONZERO COUNTS. */
/*       --------------------------------------------------------- */
    *nlnz = 0;
    i__1 = *neqns;
    for (k = 1; k <= i__1; ++k) {
	temp = colcnt[k] + weight[k];
	colcnt[k] = temp;
	*nlnz += temp;
	parent = etpar[k];
	if (parent != 0) {
	    colcnt[parent] += temp;
	}
/* L700: */
    }

    return 0;
} /* fcnthn_ */

