/* blkfct.f -- translated by f2c (version 20230428).
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
/*   Last modified:  March 6, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********     BLKFCT .....  BLOCK GENERAL SPARSE LDL'         ********* */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE CALLS THE BLOCK GENERAL SPARSE LDL' ROUTINE, */
/*       BLKFC2. */

/*   INPUT PARAMETERS: */
/*       OUTUNT          -   OUTPUT UNIT. */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       NUNROL          -   LOOP UNROLLING LEVEL. */
/*       XSUPER          -   SUPERNODE PARTITION. */
/*       SNODE           -   MAPS EACH COLUMN TO THE SUPERNODE CONTAINING */
/*                           IT. */
/*       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT */
/*                           INTO CACHE. */
/*       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE (INCLUDING */
/*                           THE DIAGONAL ELEMENTS). */
/*       (XLNZ,LNZ)      -   ON INPUT, CONTAINS MATRIX TO BE FACTORED. */
/*       IWSIZ           -   SIZE OF INTEGER WORKING STORAGE */
/*       TMPSIZ          -   SIZE OF FLOATING POINT WORKING STORAGE. */
/*       MMPYN           -   EXTERNAL ROUTINE: MATRIX-MATRIX MULTIPLY. */
/*       SMXPY           -   EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY. */

/*   OUTPUT PARAMETERS: */
/*       LNZ             -   ON OUTPUT, CONTAINS TRIANGULAR FACTOR */
/*                           AND DIAGONAL MATRIX OF LDL' DECOMPOSITION. */
/*       DIAG            -   DIAGONAL MATRIX OF LDL' DECOMPOSITION. */
/*       IFLAG           -   ERROR FLAG. */
/*                               0: SUCCESSFUL FACTORIZATION. */
/*                              -1: ZERO DIAGONAL ENCOUNTERED, */
/*                                  MATRIX IS SINGULAR. */
/*                              -2: INSUFFICIENT WORKING STORAGE */
/*                                  [TEMP(*)]. */
/*                              -3: INSUFFICIENT WORKING STORAGE */
/*                                  [IWORK(*)]. */

/*   WORKING PARAMETERS: */
/*       IWORK           -   INTEGER WORKING STORAGE OF LENGTH */
/*                           2*NEQNS + 2*NSUPER. */
/*       TMPVEC          -   DOUBLE PRECISION WORKING STORAGE OF LENGTH */
/*                           NEQNS. */

/* *********************************************************************** */

/* Subroutine */ int blkfct_(integer *neqns, integer *nsuper,
	 integer *nunrol, integer *xsuper, integer *snode, integer *split, 
	integer *xlindx, integer *lindx, integer *xlnz, doublereal *lnz, 
	doublereal *diag, integer *iwsiz, integer *iwork, integer *tmpsiz, 
	doublereal *tmpvec, integer *iflag)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer jcol;
    extern int mmpy1_(integer *, integer *, integer *, integer *,
         doublereal *, doublereal *, integer *);
    extern int mmpy2_(integer *, integer *, integer *, integer *,
         doublereal *, doublereal *, integer *);
    extern int mmpy4_(integer *, integer *, integer *, integer *,
         doublereal *, doublereal *, integer *);
    extern int mmpy8_(integer *, integer *, integer *, integer *,
         doublereal *, doublereal *, integer *);

    extern /* Subroutine */ int blkfc2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, int(integer *, integer *, integer *, integer *,
         doublereal *, doublereal *, integer *), int(integer *, integer *, doublereal *, integer *,
            doublereal *));
    extern int smxpy1_(integer *, integer *, doublereal *, integer *, 
            doublereal *);
    extern int smxpy2_(integer *, integer *, doublereal *, integer *, 
            doublereal *);
    extern int smxpy4_(integer *, integer *, doublereal *, integer *, 
            doublereal *);
    extern int smxpy8_(integer *, integer *, doublereal *, integer *, 
            doublereal *);



/* *********************************************************************** */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */


/*       --------------- */
/*       LOCAL VARIABLE. */
/*       --------------- */

/* ********************************************************************* */

    /* Parameter adjustments */
    --tmpvec;
    --iwork;
    --diag;
    --lnz;
    --xlnz;
    --lindx;
    --xlindx;
    --split;
    --snode;
    --xsuper;

    /* Function Body */
    *iflag = 0;
    if (*iwsiz < (*neqns << 1) + (*nsuper << 1)) {
	*iflag = -3;
	return 0;
    }
    *iflag = 0;
    if (*iwsiz < (*neqns << 1) + (*nsuper << 1)) {
	*iflag = -3;
	return 0;
    }
/* hao_begin */
    if (*nunrol == 8) {
	blkfc2_(nsuper, &xsuper[1], &snode[1], &split[1], &xlindx[1], &lindx[
		1], &xlnz[1], &lnz[1], &iwork[1], &iwork[*nsuper + 1], &iwork[
		(*nsuper << 1) + 1], &iwork[(*nsuper << 1) + *neqns + 1], 
		tmpsiz, &tmpvec[1], iflag, &mmpy8_, &smxpy8_);
    } else if (*nunrol == 4) {
	blkfc2_(nsuper, &xsuper[1], &snode[1], &split[1], &xlindx[1], &lindx[
		1], &xlnz[1], &lnz[1], &iwork[1], &iwork[*nsuper + 1], &iwork[
		(*nsuper << 1) + 1], &iwork[(*nsuper << 1) + *neqns + 1], 
		tmpsiz, &tmpvec[1], iflag, &mmpy4_, &smxpy4_);
    } else if (*nunrol == 2) {
	blkfc2_(nsuper, &xsuper[1], &snode[1], &split[1], &xlindx[1], &lindx[
		1], &xlnz[1], &lnz[1], &iwork[1], &iwork[*nsuper + 1], &iwork[
		(*nsuper << 1) + 1], &iwork[(*nsuper << 1) + *neqns + 1], 
		tmpsiz, &tmpvec[1], iflag, &mmpy2_, &smxpy2_);
    } else {
	blkfc2_(nsuper, &xsuper[1], &snode[1], &split[1], &xlindx[1], &lindx[
		1], &xlnz[1], &lnz[1], &iwork[1], &iwork[*nsuper + 1], &iwork[
		(*nsuper << 1) + 1], &iwork[(*nsuper << 1) + *neqns + 1], 
		tmpsiz, &tmpvec[1], iflag, &mmpy1_, &smxpy1_);
    }
/* hao_end */
    if (*iflag == -1) {
	return 0;
    } else if (*iflag == -2) {
	return 0;
    }
    if (*iflag == 0) {
	i__1 = *neqns;
	for (jcol = 1; jcol <= i__1; ++jcol) {
	    diag[jcol] = lnz[xlnz[jcol]];
/* L100: */
	}
    }
    return 0;
} /* blkfct_ */

