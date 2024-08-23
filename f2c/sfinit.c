/* sfinit.f -- translated by f2c (version 20230428).
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
#include <stdio.h>

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  January 12, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **************    SFINIT  ..... SET UP FOR SYMB. FACT.     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE COMPUTES THE STORAGE REQUIREMENTS AND SETS UP */
/*       PRELIMINARY DATA STRUCTURES FOR THE SYMBOLIC FACTORIZATION. */

/*   NOTE: */
/*       THIS VERSION PRODUCES THE MAXIMAL SUPERNODE PARTITION (I.E., */
/*       THE ONE WITH THE FEWEST POSSIBLE SUPERNODES). */

/*   INPUT PARAMETERS: */
/*       OUTUNT      -   OUTPUT UNIT. */
/*       NEQNS       -   NUMBER OF EQUATIONS. */
/*       NNZA        -   LENGTH OF ADJACENCY STRUCTURE. */
/*       XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                       TO THE ADJACENCY STRUCTURE. */
/*       ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING */
/*                       THE ADJACENCY STRUCTURE. */
/*       PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                       POSTORDERING. */
/*       INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE */
/*                       INVERSE OF THE POSTORDERING. */
/*       IWSIZ       -   SIZE OF INTEGER WORKING STORAGE. */

/*   OUTPUT PARAMETERS: */
/*       COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER */
/*                       OF NONZEROS IN EACH COLUMN OF THE FACTOR, */
/*                       INCLUDING THE DIAGONAL ENTRY. */
/*       NNZL        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING */
/*                       THE DIAGONAL ENTRIES. */
/*       NSUB        -   NUMBER OF SUBSCRIPTS. */
/*       NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS). */
/*       SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING */
/*                       SUPERNODE MEMBERSHIP. */
/*       XSUPER(*)   -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE */
/*                       SUPERNODE PARTITIONING. */
/*       IFLAG(*)    -   ERROR FLAG. */
/*                          0: SUCCESSFUL SF INITIALIZATION. */
/*                         -1: INSUFFICENT WORKING STORAGE */
/*                             [IWORK(*)]. */

/*   WORK PARAMETERS: */
/*       IWORK(*)    -   INTEGER WORK ARRAY OF LENGTH 7*NEQNS+3. */

/*   FIRST CREATED ON    NOVEMEBER 14, 1994. */
/*   LAST UPDATED ON     January 12, 1995. */

/* *********************************************************************** */

/* Subroutine */ int sfinit_(integer *neqns, integer *nnza,
	integer *xadj, integer *adjncy, integer *perm, integer *invp, integer 
	*colcnt, integer *nnzl, integer *nsub, integer *nsuper, integer *
	snode, integer *xsuper, integer *iwsiz, integer *iwork,
	int *iflag)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int fsup1_(integer *, integer *, integer *, 
	    integer *, integer *, integer *), fsup2_(integer *, integer *, 
	    integer *, integer *, integer *), fcnthn_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), chordr_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), etordr_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *);

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/* *********************************************************************** */

/*       -------------------------------------------------------- */
/*       RETURN IF THERE IS INSUFFICIENT INTEGER WORKING STORAGE. */
/*       -------------------------------------------------------- */
    /* Parameter adjustments */
    --iwork;
    --xsuper;
    --snode;
    --colcnt;
    --invp;
    --perm;
    --xadj;
    --adjncy;
    printf("after shift\n");

    /* Function Body */
    *iflag = 0;
    printf("iflag %d\n", *iflag);
    printf("iwsiz %d\n", *iwsiz);
    printf("neqns %d\n", *neqns);
    if (*iwsiz < *neqns * 7 + 3) {
	*iflag = -1;
	return 0;
    }

/*       ------------------------------------------ */
/*       COMPUTE ELIMINATION TREE AND POSTORDERING. */
/*       ------------------------------------------ */
printf("brgotr etordr_\n");
    etordr_(neqns, &xadj[1], &adjncy[1], &perm[1], &invp[1], &iwork[1], &
	    iwork[*neqns + 1], &iwork[(*neqns << 1) + 1], &iwork[*neqns * 3 + 
	    1]);

/*       --------------------------------------------- */
/*       COMPUTE ROW AND COLUMN FACTOR NONZERO COUNTS. */
/*       --------------------------------------------- */
printf("brgotr fcnthn_\n");
    fcnthn_(neqns, nnza, &xadj[1], &adjncy[1], &perm[1], &invp[1], &iwork[1], 
	    &snode[1], &colcnt[1], nnzl, &iwork[*neqns + 1], &iwork[(*neqns <<
	     1) + 1], &xsuper[1], &iwork[*neqns * 3 + 1], &iwork[(*neqns << 2)
	     + 2], &iwork[*neqns * 5 + 3], &iwork[*neqns * 6 + 4]);

/*       --------------------------------------------------------- */
/*       REARRANGE CHILDREN SO THAT THE LAST CHILD HAS THE MAXIMUM */
/*       NUMBER OF NONZEROS IN ITS COLUMN OF L. */
/*       --------------------------------------------------------- */
printf("brgotr chordr_\n");
    chordr_(neqns, &xadj[1], &adjncy[1], &perm[1], &invp[1], &colcnt[1], &
	    iwork[1], &iwork[*neqns + 1], &iwork[(*neqns << 1) + 1], &iwork[*
	    neqns * 3 + 1]);

/*       ---------------- */
/*       FIND SUPERNODES. */
/*       ---------------- */
printf("brgotr fsup1\n");
    fsup1_(neqns, &iwork[1], &colcnt[1], nsub, nsuper, &snode[1]);
printf("brgotr fsup2\n");
    fsup2_(neqns, nsuper, &iwork[1], &snode[1], &xsuper[1]);

    return 0;
} /* sfinit_ */

