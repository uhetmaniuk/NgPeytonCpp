 /* chlsup.f -- translated by f2c (version 20230428).
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
/* ******     CHLSUP .... DENSE CHOLESKY WITHIN SUPERNODE   ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*     PURPOSE - THIS ROUTINE PERFORMS CHOLESKY */
/*               FACTORIZATION ON THE COLUMNS OF A SUPERNODE */
/*               THAT HAVE RECEIVED ALL UPDATES FROM COLUMNS */
/*               EXTERNAL TO THE SUPERNODE. */

/*     INPUT PARAMETERS - */
/*        M      - NUMBER OF ROWS (LENGTH OF THE FIRST COLUMN). */
/*        N      - NUMBER OF COLUMNS IN THE SUPERNODE. */
/*        XPNT   - XPNT(J+1) POINTS ONE LOCATION BEYOND THE END */
/*                 OF THE J-TH COLUMN OF THE SUPERNODE. */
/*        X(*)   - CONTAINS THE COLUMNS OF OF THE SUPERNODE TO */
/*                 BE FACTORED. */
/*        SMXPY  - EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY. */

/*     OUTPUT PARAMETERS - */
/*        X(*)   - ON OUTPUT, CONTAINS THE FACTORED COLUMNS OF */
/*                 THE SUPERNODE. */
/*        IFLAG  - UNCHANGED IF THERE IS NO ERROR. */
/*                 =1 IF NONPOSITIVE DIAGONAL ENTRY IS ENCOUNTERED. */

/* *********************************************************************** */

/* Subroutine */ int chlsup_(integer *m, integer *n, integer *split, integer *
	xpnt, doublereal *x, integer *iflag, 
int (*mmpyn)(integer *, integer *, integer *, integer *,
         doublereal *, doublereal *, integer *), int (*smxpy)(integer *, integer *, doublereal *, integer *,
            doublereal *))
{
    static integer q, mm, nn, jblk, jpnt;
    extern /* Subroutine */ int pchol_(integer *, integer *, integer *, 
	    doublereal *, integer *, int(integer *, integer *, doublereal *, integer *, 
            doublereal *));
    static integer fstcol, nxtcol;


/* *********************************************************************** */

/*     ----------- */
/*     PARAMETERS. */
/*     ----------- */





/*     ---------------- */
/*     LOCAL VARIABLES. */
/*     ---------------- */


/* *********************************************************************** */

    /* Parameter adjustments */
    --x;
    --xpnt;
    --split;

    /* Function Body */
    jblk = 0;
    fstcol = 1;
    mm = *m;
    jpnt = xpnt[fstcol];

/*       ---------------------------------------- */
/*       FOR EACH BLOCK JBLK IN THE SUPERNODE ... */
/*       ---------------------------------------- */
L100:
    if (fstcol <= *n) {
	++jblk;
	nn = split[jblk];
/*           ------------------------------------------ */
/*           ... PERFORM PARTIAL CHOLESKY FACTORIZATION */
/*               ON THE BLOCK. */
/*           ------------------------------------------ */
	pchol_(&mm, &nn, &xpnt[fstcol], &x[1], iflag, smxpy);
	if (*iflag == 1) {
	    return 0;
	}
/*           ---------------------------------------------- */
/*           ... APPLY THE COLUMNS IN JBLK TO ANY COLUMNS */
/*               OF THE SUPERNODE REMAINING TO BE COMPUTED. */
/*           ---------------------------------------------- */
	nxtcol = fstcol + nn;
	q = *n - nxtcol + 1;
	mm -= nn;
	jpnt = xpnt[nxtcol];
	if (q > 0) {
	    (*mmpyn)(&mm, &nn, &q, &xpnt[fstcol], &x[1], &x[jpnt], &mm);
	}
	fstcol = nxtcol;
	goto L100;
    }

    return 0;
} /* chlsup_ */

