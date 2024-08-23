/* genmmd.f -- translated by f2c (version 20230428).
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
/* --- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = GENMMD */
/*  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984 */
/* *********************************************************************** */
/* *********************************************************************** */
/* ****     GENMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE IMPLEMENTS THE MINIMUM DEGREE */
/*        ALGORITHM.  IT MAKES USE OF THE IMPLICIT REPRESENTATION */
/*        OF ELIMINATION GRAPHS BY QUOTIENT GRAPHS, AND THE */
/*        NOTION OF INDISTINGUISHABLE NODES.  IT ALSO IMPLEMENTS */
/*        THE MODIFICATIONS BY MULTIPLE ELIMINATION AND MINIMUM */
/*        EXTERNAL DEGREE. */
/*        --------------------------------------------- */
/*        CAUTION - THE ADJACENCY VECTOR ADJNCY WILL BE */
/*        DESTROYED. */
/*        --------------------------------------------- */

/*     INPUT PARAMETERS - */
/*        NEQNS  - NUMBER OF EQUATIONS. */
/*        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE. */
/*        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION. */
/*        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER */
/*                 (ANY SMALLER ESTIMATE WILL DO) FOR MARKING */
/*                 NODES. */

/*     OUTPUT PARAMETERS - */
/*        PERM   - THE MINIMUM DEGREE ORDERING. */
/*        INVP   - THE INVERSE OF PERM. */
/*        NNZL   - NUMBER OF NONZEROS IN LOWER TRIAGULAR FACTOR. */
/*        NOFSUB - NUMBER OF SUBSCRIPTS FOR THE COMPRESSED STORAGE */
/*                 SCHEME. */
/*        COLCNT - NUMBER OF NONZEROS IN EACH FACTOR COLUMN, INCLUDING */
/*                 THE DIAGONAL ENTRY (DEGREE+1). */
/*        NSUPER - NUMBER OF SUPERNODES. */
/*        XSUPER - FIRST COLUMN OF EACH SUPERNODE. */
/*        SNODE  - SUPERNODE MEMBERSHIP OF EACH COLUMN. */

/*     WORKING PARAMETERS - */
/*        DHEAD  - VECTOR FOR HEAD OF DEGREE LISTS. */
/*        INVP   - USED TEMPORARILY FOR DEGREE FORWARD LINK. */
/*        PERM   - USED TEMPORARILY FOR DEGREE BACKWARD LINK. */
/*        QSIZE  - VECTOR FOR SIZE OF SUPERNODES. */
/*        LLIST  - VECTOR FOR TEMPORARY LINKED LISTS. */
/*        MARKER - A TEMPORARY MARKER VECTOR. */

/*     PROGRAM SUBROUTINES - */
/*        MMDELM, MMDINT, MMDNUM, MMDUPD. */

/* *********************************************************************** */

/* Subroutine */ int genmmd_(integer *neqns, integer *xadj, integer *adjncy, 
	integer *invp, integer *perm, integer *delta, integer *dhead, integer 
	*qsize, integer *llist, integer *marker, integer *maxint, integer *
	nnzl, integer *nofsub, integer *colcnt, integer *nsuper, integer *
	xsuper, integer *snode)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, cc, tag, num, mdeg, kcol, ehead, mdlmt, mdnode;
    extern /* Subroutine */ int mmdelm_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), mmdupd_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), mmdint_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);
    static integer fstcol;
    extern /* Subroutine */ int mmdnum_(integer *, integer *, integer *, 
	    integer *);
    static integer nextmd, lstcol, ksuper;


/* *********************************************************************** */


/* *********************************************************************** */

    /* Parameter adjustments */
    --snode;
    --xsuper;
    --colcnt;
    --marker;
    --llist;
    --qsize;
    --dhead;
    --perm;
    --invp;
    --adjncy;
    --xadj;

    /* Function Body */
    if (*neqns <= 0) {
	return 0;
    }

/*        ------------------------------------------------ */
/*        INITIALIZATION FOR THE MINIMUM DEGREE ALGORITHM. */
/*        ------------------------------------------------ */
    *nsuper = 0;
    mmdint_(neqns, &xadj[1], &adjncy[1], &dhead[1], &invp[1], &perm[1], &
	    qsize[1], &llist[1], &marker[1]);

/*        ---------------------------------------------- */
/*        NUM COUNTS THE NUMBER OF ORDERED NODES PLUS 1. */
/*        ---------------------------------------------- */
    num = 1;

/*        ----------------------------- */
/*        ELIMINATE ALL ISOLATED NODES. */
/*        ----------------------------- */
    nextmd = dhead[1];
L100:
    if (nextmd <= 0) {
	goto L200;
    }
    mdnode = nextmd;
    nextmd = invp[mdnode];
    marker[mdnode] = *maxint;
    invp[mdnode] = -num;
    ++(*nsuper);
    xsuper[*nsuper] = num;
    colcnt[num] = 1;
    ++num;
    goto L100;

L200:
/*        ---------------------------------------- */
/*        SEARCH FOR NODE OF THE MINIMUM DEGREE. */
/*        MDEG IS THE CURRENT MINIMUM DEGREE; */
/*        TAG IS USED TO FACILITATE MARKING NODES. */
/*        ---------------------------------------- */
    if (num > *neqns) {
	goto L1000;
    }
    tag = 1;
    dhead[1] = 0;
    mdeg = 2;
L300:
    if (dhead[mdeg] > 0) {
	goto L400;
    }
    ++mdeg;
    goto L300;
L400:
/*            ------------------------------------------------- */
/*            USE VALUE OF DELTA TO SET UP MDLMT, WHICH GOVERNS */
/*            WHEN A DEGREE UPDATE IS TO BE PERFORMED. */
/*            ------------------------------------------------- */
    mdlmt = mdeg + *delta;
    ehead = 0;

L500:
    mdnode = dhead[mdeg];
    if (mdnode > 0) {
	goto L600;
    }
    ++mdeg;
    if (mdeg > mdlmt) {
	goto L900;
    }
    goto L500;
L600:
/*                ---------------------------------------- */
/*                REMOVE MDNODE FROM THE DEGREE STRUCTURE. */
/*                ---------------------------------------- */
    nextmd = invp[mdnode];
    dhead[mdeg] = nextmd;
    if (nextmd > 0) {
	perm[nextmd] = -mdeg;
    }
    ++(*nsuper);
    xsuper[*nsuper] = num;
    colcnt[num] = mdeg + qsize[mdnode] - 1;
    invp[mdnode] = -num;
    if (num + qsize[mdnode] > *neqns) {
	goto L1000;
    }
/*                ---------------------------------------------- */
/*                ELIMINATE MDNODE AND PERFORM QUOTIENT GRAPH */
/*                TRANSFORMATION.  RESET TAG VALUE IF NECESSARY. */
/*                ---------------------------------------------- */
    ++tag;
    if (tag < *maxint) {
	goto L800;
    }
    tag = 1;
    i__1 = *neqns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (marker[i__] < *maxint) {
	    marker[i__] = 0;
	}
/* L700: */
    }
L800:
    mmdelm_(&mdnode, &xadj[1], &adjncy[1], &dhead[1], &invp[1], &perm[1], &
	    qsize[1], &llist[1], &marker[1], maxint, &tag);
    num += qsize[mdnode];
    llist[mdnode] = ehead;
    ehead = mdnode;
    if (*delta >= 0) {
	goto L500;
    }
L900:
/*            ------------------------------------------- */
/*            UPDATE DEGREES OF THE NODES INVOLVED IN THE */
/*            MINIMUM DEGREE NODES ELIMINATION. */
/*            ------------------------------------------- */
    if (num > *neqns) {
	goto L1000;
    }
    mmdupd_(&ehead, neqns, &xadj[1], &adjncy[1], delta, &mdeg, &dhead[1], &
	    invp[1], &perm[1], &qsize[1], &llist[1], &marker[1], maxint, &tag)
	    ;
    goto L300;

L1000:
    xsuper[*nsuper + 1] = *neqns + 1;
    *nnzl = 0;
    *nofsub = 0;
    i__1 = *nsuper;
    for (ksuper = 1; ksuper <= i__1; ++ksuper) {
	fstcol = xsuper[ksuper];
	lstcol = xsuper[ksuper + 1] - 1;
	cc = colcnt[fstcol];
	*nofsub += cc;
	i__2 = lstcol;
	for (kcol = fstcol; kcol <= i__2; ++kcol) {
	    snode[kcol] = ksuper;
	    colcnt[kcol] = cc;
	    *nnzl += cc;
	    --cc;
/* L1100: */
	}
/* L1200: */
    }
    mmdnum_(neqns, &perm[1], &invp[1], &qsize[1]);
    return 0;

} /* genmmd_ */

