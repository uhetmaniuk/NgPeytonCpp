/* blkfc2.f -- translated by f2c (version 20230428).
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
/*   Last modified:  March 6, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *********     BLKFC2 .....  BLOCK GENERAL SPARSE CHOLESKY     ********* */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE FACTORS A SPARSE POSITIVE DEFINITE MATRIX. */
/*       THE COMPUTATION IS ORGANIZED AROUND KERNELS THAT PERFORM */
/*       SUPERNODE-TO-SUPERNODE UPDATES, I.E., BLOCK-TO-BLOCK UPDATES. */

/*   INPUT PARAMETERS: */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   SUPERNODE PARTITION. */
/*       SNODE           -   MAPS EACH COLUMN TO THE SUPERNODE CONTAINING */
/*                           IT. */
/*       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT */
/*                           INTO CACHE. */
/*       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE (INCLUDING */
/*                           THE DIAGONAL ELEMENTS). */
/*       (XLNZ,LNZ)      -   ON INPUT, CONTAINS MATRIX TO BE FACTORED. */
/*       TMPSIZ          -   SIZE OF TEMPORARY WORKING STORAGE. */
/*       MMPYN           -   EXTERNAL ROUTINE: MATRIX-MATRIX MULTIPLY. */
/*       SMXPY           -   EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY. */

/*   OUTPUT PARAMETERS: */
/*       LNZ             -   ON OUTPUT, CONTAINS CHOLESKY FACTOR. */
/*       IFLAG           -   ERROR FLAG. */
/*                               0: SUCCESSFUL FACTORIZATION. */
/*                              -1: NONPOSITIVE DIAGONAL ENCOUNTERED, */
/*                                  MATRIX IS NOT POSITIVE DEFINITE. */
/*                              -2: INSUFFICIENT WORKING STORAGE */
/*                                  [TEMP(*)]. */

/*   WORKING PARAMETERS: */
/*       LINK            -   LINKS TOGETHER THE SUPERNODES IN A SUPERNODE */
/*                           ROW. */
/*       LENGTH          -   LENGTH OF THE ACTIVE PORTION OF EACH */
/*                           SUPERNODE. */
/*       INDMAP          -   VECTOR OF SIZE NEQNS INTO WHICH THE GLOBAL */
/*                           INDICES ARE SCATTERED. */
/*       RELIND          -   MAPS LOCATIONS IN THE UPDATING COLUMNS TO */
/*                           THE CORRESPONDING LOCATIONS IN THE UPDATED */
/*                           COLUMNS.  (RELIND IS GATHERED FROM INDMAP). */
/*       TEMP            -   REAL VECTOR FOR ACCUMULATING UPDATES.  MUST */
/*                           ACCOMODATE ALL COLUMNS OF A SUPERNODE. */

/* *********************************************************************** */

/* Subroutine */ int blkfc2_(integer *nsuper, integer *xsuper, integer *snode,
	 integer *split, integer *xlindx, integer *lindx, integer *xlnz, 
	doublereal *lnz, integer *link, integer *length, integer *indmap, 
	integer *relind, integer *tmpsiz, doublereal *temp, integer *iflag, 
	int (*mmpyn)(integer *, integer *, integer *, integer *,
         doublereal *, doublereal *, integer *), int (*smxpy)(integer *, integer *, doublereal *, integer *,
            doublereal *))
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, ilen, jlen, klen, jsup, ksup;
    extern /* Subroutine */ int mmpy_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
int(integer *, integer *, integer *, integer *,
         doublereal *, doublereal *, integer *))
	    ;
    static integer fjcol, fkcol, ljcol;
    extern /* Subroutine */ int assmb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);
    static integer klast, kdpnt, ilpnt, jlpnt, klpnt, store;
    extern /* Subroutine */ int mmpyi_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    static integer jxpnt, kxpnt, inddif;
    extern /* Subroutine */ int igathr_(integer *, integer *, integer *, 
	    integer *), ldindx_(integer *, integer *, integer *);
    static integer njcols, nkcols;
    extern /* Subroutine */ int chlsup_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, int(integer *, integer *, integer *, integer *,
         doublereal *, doublereal *, integer *), int(integer *, integer *, doublereal *, integer *,        
            doublereal *));
    static integer ncolup, kfirst, nxtcol, nxksup, nxtsup;



/* ********************************************************************* */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */


/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */


/* ********************************************************************* */

    /* Parameter adjustments */
    --temp;
    --relind;
    --indmap;
    --length;
    --link;
    --lnz;
    --xlnz;
    --lindx;
    --xlindx;
    --split;
    --snode;
    --xsuper;

    /* Function Body */
    *iflag = 0;

/*       ----------------------------------------------------------- */
/*       INITIALIZE EMPTY ROW LISTS IN LINK(*) AND ZERO OUT TEMP(*). */
/*       ----------------------------------------------------------- */
    i__1 = *nsuper;
    for (jsup = 1; jsup <= i__1; ++jsup) {
	link[jsup] = 0;
/* L100: */
    }
    i__1 = *tmpsiz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp[i__] = 0.;
/* L200: */
    }

/*       --------------------------- */
/*       FOR EACH SUPERNODE JSUP ... */
/*       --------------------------- */
    i__1 = *nsuper;
    for (jsup = 1; jsup <= i__1; ++jsup) {

/*           ------------------------------------------------ */
/*           FJCOL  ...  FIRST COLUMN OF SUPERNODE JSUP. */
/*           LJCOL  ...  LAST COLUMN OF SUPERNODE JSUP. */
/*           NJCOLS ...  NUMBER OF COLUMNS IN SUPERNODE JSUP. */
/*           JLEN   ...  LENGTH OF COLUMN FJCOL. */
/*           JXPNT  ...  POINTER TO INDEX OF FIRST */
/*                       NONZERO IN COLUMN FJCOL. */
/*           ------------------------------------------------ */
	fjcol = xsuper[jsup];
	njcols = xsuper[jsup + 1] - fjcol;
	ljcol = fjcol + njcols - 1;
	jlen = xlnz[fjcol + 1] - xlnz[fjcol];
	jxpnt = xlindx[jsup];

/*           ----------------------------------------------------- */
/*           SET UP INDMAP(*) TO MAP THE ENTRIES IN UPDATE COLUMNS */
/*           TO THEIR CORRESPONDING POSITIONS IN UPDATED COLUMNS, */
/*           RELATIVE THE THE BOTTOM OF EACH UPDATED COLUMN. */
/*           ----------------------------------------------------- */
	ldindx_(&jlen, &lindx[jxpnt], &indmap[1]);

/*           ----------------------------------------- */
/*           FOR EVERY SUPERNODE KSUP IN ROW(JSUP) ... */
/*           ----------------------------------------- */
	ksup = link[jsup];
L300:
	if (ksup > 0) {
	    nxksup = link[ksup];

/*               ------------------------------------------------------- */
/*               GET INFO ABOUT THE CMOD(JSUP,KSUP) UPDATE. */

/*               FKCOL  ...  FIRST COLUMN OF SUPERNODE KSUP. */
/*               NKCOLS ...  NUMBER OF COLUMNS IN SUPERNODE KSUP. */
/*               KLEN   ...  LENGTH OF ACTIVE PORTION OF COLUMN FKCOL. */
/*               KXPNT  ...  POINTER TO INDEX OF FIRST NONZERO IN ACTIVE */
/*                           PORTION OF COLUMN FJCOL. */
/*               ------------------------------------------------------- */
	    fkcol = xsuper[ksup];
	    nkcols = xsuper[ksup + 1] - fkcol;
	    klen = length[ksup];
	    kxpnt = xlindx[ksup + 1] - klen;

/*               ------------------------------------------- */
/*               PERFORM CMOD(JSUP,KSUP), WITH SPECIAL CASES */
/*               HANDLED DIFFERENTLY. */
/*               ------------------------------------------- */

	    if (klen != jlen) {

/*                   ------------------------------------------- */
/*                   SPARSE CMOD(JSUP,KSUP). */

/*                   NCOLUP ... NUMBER OF COLUMNS TO BE UPDATED. */
/*                   ------------------------------------------- */

		i__2 = klen - 1;
		for (i__ = 0; i__ <= i__2; ++i__) {
		    nxtcol = lindx[kxpnt + i__];
		    if (nxtcol > ljcol) {
			goto L500;
		    }
/* L400: */
		}
		i__ = klen;
L500:
		ncolup = i__;

		if (nkcols == 1) {

/*                       ---------------------------------------------- */
/*                       UPDATING TARGET SUPERNODE BY TRIVIAL */
/*                       SUPERNODE (WITH ONE COLUMN). */

/*                       KLPNT  ...  POINTER TO FIRST NONZERO IN ACTIVE */
/*                                   PORTION OF COLUMN FKCOL. */
/*                       KDPNT  ...  POINTER TO DIAGONAL ENTRY OF */
/*                                   COLUMN FKCOL. */
/*                       ---------------------------------------------- */
		    klpnt = xlnz[fkcol + 1] - klen;
		    kdpnt = xlnz[fkcol];
		    mmpyi_(&klen, &ncolup, &lindx[kxpnt], &lnz[klpnt], &lnz[
			    kdpnt], &xlnz[1], &lnz[1], &indmap[1]);

		} else {

/*                       -------------------------------------------- */
/*                       KFIRST ...  FIRST INDEX OF ACTIVE PORTION OF */
/*                                   SUPERNODE KSUP (FIRST COLUMN TO */
/*                                   BE UPDATED). */
/*                       KLAST  ...  LAST INDEX OF ACTIVE PORTION OF */
/*                                   SUPERNODE KSUP. */
/*                       -------------------------------------------- */

		    kfirst = lindx[kxpnt];
		    klast = lindx[kxpnt + klen - 1];
		    inddif = indmap[kfirst] - indmap[klast];

		    if (inddif < klen) {

/*                           --------------------------------------- */
/*                           DENSE CMOD(JSUP,KSUP). */

/*                           ILPNT  ...  POINTER TO FIRST NONZERO IN */
/*                                       COLUMN KFIRST. */
/*                           ILEN   ...  LENGTH OF COLUMN KFIRST. */
/*                           --------------------------------------- */
			ilpnt = xlnz[kfirst];
			ilen = xlnz[kfirst + 1] - ilpnt;
			mmpy_(&klen, &nkcols, &ncolup, &split[fkcol], &xlnz[
				fkcol], &lnz[1], &lnz[ilpnt], &ilen,
				mmpyn);

		    } else {

/*                           ------------------------------- */
/*                           GENERAL SPARSE CMOD(JSUP,KSUP). */
/*                           COMPUTE CMOD(JSUP,KSUP) UPDATE */
/*                           IN WORK STORAGE. */
/*                           ------------------------------- */
			store = klen * ncolup - ncolup * (ncolup - 1) / 2;
			if (store > *tmpsiz) {
			    *iflag = -2;
			    return 0;
			}
			mmpy_(&klen, &nkcols, &ncolup, &split[fkcol], &xlnz[
				fkcol], &lnz[1], &temp[1], &klen, mmpyn)
				;
/*                           ---------------------------------------- */
/*                           GATHER INDICES OF KSUP RELATIVE TO JSUP. */
/*                           ---------------------------------------- */
			igathr_(&klen, &lindx[kxpnt], &indmap[1], &relind[1]);
/*                           -------------------------------------- */
/*                           INCORPORATE THE CMOD(JSUP,KSUP) BLOCK */
/*                           UPDATE INTO THE TO APPROPRIATE COLUMNS */
/*                           OF L. */
/*                           -------------------------------------- */
			assmb_(&klen, &ncolup, &temp[1], &relind[1], &xlnz[
				fjcol], &lnz[1], &jlen);

		    }

		}

	    } else {

/*                   ---------------------------------------------- */
/*                   DENSE CMOD(JSUP,KSUP). */
/*                   JSUP AND KSUP HAVE IDENTICAL STRUCTURE. */

/*                   JLPNT  ...  POINTER TO FIRST NONZERO IN COLUMN */
/*                               FJCOL. */
/*                   ---------------------------------------------- */
		jlpnt = xlnz[fjcol];
		mmpy_(&klen, &nkcols, &njcols, &split[fkcol], &xlnz[fkcol], &
			lnz[1], &lnz[jlpnt], &jlen, mmpyn);
		ncolup = njcols;
		if (klen > njcols) {
		    nxtcol = lindx[jxpnt + njcols];
		}

	    }

/*               ------------------------------------------------ */
/*               LINK KSUP INTO LINKED LIST OF THE NEXT SUPERNODE */
/*               IT WILL UPDATE AND DECREMENT KSUP'S ACTIVE */
/*               LENGTH. */
/*               ------------------------------------------------ */
	    if (klen > ncolup) {
		nxtsup = snode[nxtcol];
		link[ksup] = link[nxtsup];
		link[nxtsup] = ksup;
		length[ksup] = klen - ncolup;
	    } else {
		length[ksup] = 0;
	    }

/*               ------------------------------- */
/*               NEXT UPDATING SUPERNODE (KSUP). */
/*               ------------------------------- */
	    ksup = nxksup;
	    goto L300;

	}

/*           ---------------------------------------------- */
/*           APPLY PARTIAL CHOLESKY TO THE COLUMNS OF JSUP. */
/*           ---------------------------------------------- */
	chlsup_(&jlen, &njcols, &split[fjcol], &xlnz[fjcol], &lnz[1], iflag,
		mmpyn, smxpy);
	if (*iflag != 0) {
	    *iflag = -1;
	    return 0;
	}

/*           ----------------------------------------------- */
/*           INSERT JSUP INTO LINKED LIST OF FIRST SUPERNODE */
/*           IT WILL UPDATE. */
/*           ----------------------------------------------- */
	if (jlen > njcols) {
	    nxtcol = lindx[jxpnt + njcols];
	    nxtsup = snode[nxtcol];
	    link[jsup] = link[nxtsup];
	    link[nxtsup] = jsup;
	    length[jsup] = jlen - njcols;
	} else {
	    length[jsup] = 0;
	}

/* L600: */
    }

    return 0;
} /* blkfc2_ */

