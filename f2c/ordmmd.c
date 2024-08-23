/* ordmmd.f -- translated by f2c (version 20230428).
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
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****     ORDMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE CALLS LIU'S MULTIPLE MINIMUM DEGREE */
/*               ROUTINE. */

/*     INPUT PARAMETERS - */
/*        OUTUNT - OUTPUT UNIT. */
/*        NEQNS  - NUMBER OF EQUATIONS. */
/*        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE. */
/*        IWSIZ  - SIZE OF INTEGER WORKING STORAGE. */

/*     OUTPUT PARAMETERS - */
/*        PERM   - THE MINIMUM DEGREE ORDERING. */
/*        INVP   - THE INVERSE OF PERM. */
/*        NNZL   - NUMBER OF NONZEROS IN LOWER TRIAGULAR FACTOR. */
/*        NOFSUB - NUMBER OF SUBSCRIPTS FOR THE COMPRESSED STORAGE */
/*                 SCHEME. */
/*        NOFSUB - NUMBER OF SUBSCRIPTS FOR THE COMPRESSED STORAGE */
/*                 SCHEME. */
/*        COLCNT - NUMBER OF NONZEROS IN EACH FACTOR COLUMN, INCLUDING */
/*                 THE DIAGONAL ENTRY (DEGREE+1). */
/*        NSUPER - NUMBER OF SUPERNODES. */
/*        XSUPER - FIRST COLUMN OF EACH SUPERNODE. */
/*        SNODE  - SUPERNODE MEMBERSHIP OF EACH COLUMN. */
/*        SFIFLG - SFIFLG=.F. MEANS SKIP SYMBOLIC FACTORIZATION */
/*                 INITIALIZATION (SFINIT), SFIFLG=.T. MEANS EXECUTE */
/*                 SFINIT. */
/*        IFLAG  - ERROR FLAG. */
/*                   0: SUCCESSFUL ORDERING */
/*                  -1: INSUFFICIENT WORKING STORAGE */
/*                      [IWORK(*)]. */

/*     WORKING PARAMETERS - */
/*        IWORK  - INTEGER WORKSPACE OF LENGTH 4*NEQNS. */

/* *********************************************************************** */

/* Subroutine */ int ordmmd_(integer *neqns, integer *xadj, 
	integer *adjncy, integer *invp, integer *perm, integer *iwsiz, 
	integer *iwork, integer *nnzl, integer *nofsub, integer *colcnt, 
	integer *nsuper, integer *xsuper, integer *snode, logical *sfiflg, 
	integer *iflag)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer delta;
    extern /* Subroutine */ int genmmd_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    static integer maxint;

/* *********************************************************************** */


/* ********************************************************************* */

    /* Parameter adjustments */
    --snode;
    --xsuper;
    --colcnt;
    --iwork;
    --perm;
    --invp;
    --adjncy;
    --xadj;

    /* Function Body */
    *iflag = 0;
    if (*iwsiz < *neqns << 2) {
	*iflag = -1;
	return 0;
    }

/*       DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION. */
/*       MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER */
/*                (ANY SMALLER ESTIMATE WILL DO) FOR MARKING */
/*                NODES. */

    delta = 0;
    maxint = 32767;
    genmmd_(neqns, &xadj[1], &adjncy[1], &invp[1], &perm[1], &delta, &iwork[1]
	    , &iwork[*neqns + 1], &iwork[(*neqns << 1) + 1], &iwork[*neqns * 
	    3 + 1], &maxint, nnzl, nofsub, &colcnt[1], nsuper, &xsuper[1], &
	    snode[1]);
    *sfiflg = FALSE_;
    return 0;

} /* ordmmd_ */

