/* symfct.f -- translated by f2c (version 20230428).
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

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  February 13, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     SYMFCT ..... SYMBOLIC FACTORIZATION    ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS ROUTINE CALLS SYMFC2 WHICH PERFORMS SUPERNODAL SYMBOLIC */
/*       FACTORIZATION ON A REORDERED LINEAR SYSTEM. */

/*   INPUT PARAMETERS: */
/*       (I) OUTUNT      -   OUTPUT UNIT. */
/*       (I) NEQNS       -   NUMBER OF EQUATIONS */
/*       (I) ADJLEN      -   LENGTH OF THE ADJACENCY LIST. */
/*       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1 CONTAINING POINTERS */
/*                           TO THE ADJACENCY STRUCTURE. */
/*       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1 CONTAINING */
/*                           THE ADJACENCY STRUCTURE. */
/*       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           POSTORDERING. */
/*       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE */
/*                           INVERSE OF THE POSTORDERING. */
/*       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                           OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                           INCLUDING THE DIAGONAL ENTRY. */
/*       (I) NSUPER      -   NUMBER OF SUPERNODES. */
/*       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE */
/*                           FIRST COLUMN OF EACH SUPERNODE. */
/*       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                           SUPERNODE MEMBERSHIP. */
/*       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO BE STORED IN */
/*                           LINDX(*). */
/*       (I) IWSIZ       -   SIZE OF INTEGER WORKING STORAGE. */

/*   OUTPUT PARAMETERS: */
/*       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           INTO THE SUBSCRIPT VECTOR. */
/*       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE */
/*                           COMPRESSED SUBSCRIPTS. */
/*       (I) XLNZ        -   COLUMN POINTERS FOR L. */
/*       (I) FLAG        -   ERROR FLAG: */
/*                               0 - NO ERROR. */
/*                              -1 - INSUFFICIENT INTEGER WORKING SPACE. */
/*                              -2 - INCONSISTANCY IN THE INPUT. */

/*   WORKING PARAMETERS: */
/*       (I) IWORK       -   WORKING ARRAY OF LENGTH NSUPER+2*NEQNS. */

/* *********************************************************************** */

/* Subroutine */ int symfct_(integer *neqns, integer *adjlen,
	 integer *xadj, integer *adjncy, integer *perm, integer *invp, 
	integer *colcnt, integer *nsuper, integer *xsuper, integer *snode, 
	integer *nofsub, integer *xlindx, integer *lindx, integer *xlnz, 
	integer *iwsiz, integer *iwork, integer *flag__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int symfc2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);


/* *********************************************************************** */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/* *********************************************************************** */

    /* Parameter adjustments */
    --xlnz;
    --snode;
    --colcnt;
    --invp;
    --perm;
    --xadj;
    --adjncy;
    --iwork;
    --xlindx;
    --xsuper;
    --lindx;

    /* Function Body */
    *flag__ = 0;
    if (*iwsiz < *nsuper + (*neqns << 1) + 1) {
	*flag__ = -1;
	return 0;
    }
    symfc2_(neqns, adjlen, &xadj[1], &adjncy[1], &perm[1], &invp[1], &colcnt[
	    1], nsuper, &xsuper[1], &snode[1], nofsub, &xlindx[1], &lindx[1], 
	    &xlnz[1], &iwork[1], &iwork[*nsuper + 1], &iwork[*nsuper + *neqns 
	    + 2], flag__);
    if (*flag__ == -2) {
	return 0;
    }
    return 0;
} /* symfct_ */

