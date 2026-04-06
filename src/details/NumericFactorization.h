#ifndef NGPEYTONCPP_DETAILS_NUMERICFACTORIZATION_H
#define NGPEYTONCPP_DETAILS_NUMERICFACTORIZATION_H

#include <algorithm>
#include <iostream>

#include "details/Utilities.h"
#include "details/BLAS.h"

namespace NgPeytonCpp { namespace details { namespace f2c {

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */

/*     ------------------------------------------------------ */
/*     INPUT NUMERICAL VALUES INTO SPARSE DATA STRUCTURES ... */
/*     ------------------------------------------------------ */

template <typename Scalar, typename Index>
void inpnv(
  Index neqns, Index* __restrict xadjf, Index* __restrict adjf,
  Scalar* __restrict anzf, Index* __restrict perm, Index* __restrict invp,
  Index nsuper, Index* __restrict xsuper, Index* __restrict xlindx,
  Index* __restrict lindx, Index* __restrict xlnz, Scalar* __restrict lnz,
  Index iwsiz, Index* __restrict offset, Index& iflag) {
  /* System generated locals */
  Index i__4;

  /* Local variables */
  Index i, j, ii, jlen, oldj, last, jsuper;

  /* Function Body */
  iflag = 0;
  if (iwsiz < neqns) {
    iflag = -1;
    std::cerr << "\n";
    std::cerr << "*** INTEGER WORK SPACE = " << iwsiz << "\n";
    std::cerr << "*** IS SMALLER THAN REQUIRED = " << neqns << "\n";
    std::cerr << "\n";
    return;
  }

  /* Parameter adjustments */
  --offset;
  --lnz;
  --xlnz;
  --lindx;
  --xlindx;
  --xsuper;
  --invp;
  --perm;
  --anzf;
  --adjf;
  --xadjf;

  for (jsuper = 1; jsuper <= nsuper; ++jsuper) {
    /*           ---------------------------------------- */
    /*           FOR EACH SUPERNODE, DO THE FOLLOWING ... */
    /*           ---------------------------------------- */

    /*           ----------------------------------------------- */
    /*           FIRST GET OFFSET TO FACILITATE NUMERICAL INPUT. */
    /*           ----------------------------------------------- */
    jlen = xlindx[jsuper + 1] - xlindx[jsuper];
    for (ii = xlindx[jsuper]; ii < xlindx[jsuper + 1]; ++ii) {
      i = lindx[ii];
      --jlen;
      offset[i] = jlen;
    }

    for (j = xsuper[jsuper]; j < xsuper[jsuper + 1]; ++j) {
      /*               ----------------------------------------- */
      /*               FOR EACH COLUMN IN THE CURRENT SUPERNODE, */
      /*               FIRST INITIALIZE THE DATA STRUCTURE. */
      /*               ----------------------------------------- */
      for (ii = xlnz[j]; ii < xlnz[j + 1]; ++ii) {
        lnz[ii] = 0;
      }

      /*               ----------------------------------- */
      /*               NEXT INPUT THE INDIVIDUAL NONZEROS. */
      /*               ----------------------------------- */
      oldj = perm[j];
      last = xlnz[j + 1] - 1;
      for (ii = xadjf[oldj]; ii < xadjf[oldj + 1]; ++ii) {
        i = invp[adjf[ii]];
        if (i >= j) {
          i__4 = last - offset[i];
          lnz[i__4] = anzf[ii];
        }
      }
    }
  }
} /* inpnv */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     PCHOL .... DENSE PARTIAL CHOLESKY             ************** */
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

template <typename Scalar, typename Index>
void pchol(Index m, Index n, Index* __restrict xpnt, Scalar* __restrict x,
           Index& iflag) {
  /* Local variables */
  Index mm;
  Scalar diag;
  Index jcol, jpnt;

  /*       ------------------------------------------ */
  /*       FOR EVERY COLUMN JCOL IN THE SUPERNODE ... */
  /*       ------------------------------------------ */
  /* Parameter adjustments */
  --x;
  --xpnt;

  /* Function Body */
  mm = m;
  jpnt = xpnt[1];
  for (jcol = 1; jcol <= n; ++jcol) {
    /*           ---------------------------------- */
    /*           UPDATE JCOL WITH PREVIOUS COLUMNS. */
    /*           ---------------------------------- */
    if (jcol > 1) {
      smxpy1(mm, jcol - 1, &x[jpnt], &xpnt[1], &x[1]);
    }

    /*           --------------------------- */
    /*           COMPUTE THE DIAGONAL ENTRY. */
    /*           --------------------------- */
    diag = x[jpnt];
    if (diag == static_cast<Scalar>(0.0)) {
      iflag = 1;
      return;
    }
    diag = static_cast<Scalar>(1.0) / diag;

    /*           ---------------------------------------------------- */
    /*           SCALE COLUMN JCOL WITH RECIPROCAL OF DIAGONAL ENTRY. */
    /*           ---------------------------------------------------- */
    --mm;
    ++jpnt;
    scal(mm, diag, &x[jpnt]);
    jpnt += mm;
  }
} /* pchol */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     MMPY1  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*       LOOP UNROLLING: LEVEL 1 */

/*   INPUT PARAMETERS - */
/*       M               -   NUMBER OF ROWS IN X AND IN Y. */
/*       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS */
/*                           IN A. */
/*       Q               -   NUMBER OF COLUMNS IN A AND Y. */
/*       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE */
/*                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO */
/*                           USED TO ACCESS THE ROWS OF A. */
/*       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A. */
/*       LDY             -   LENGTH OF FIRST COLUMN OF Y. */

/*   UPDATED PARAMETERS - */
/*       Y(*)            -   ON OUTPUT, Y = Y + AX. */

/* *********************************************************************** */

template <typename Scalar, typename Index>
int mmpy1(
  Index m, Index n, Index q, Index* __restrict xpnt, Scalar* __restrict x,
  Scalar* __restrict y, Index ldy) {
  /* System generated locals */
  Index i__4, i__5, i__6;
  Scalar z__1, z__2;

  /* Local variables */
  Scalar a1;
  Index i1, j1, mm, iy, xcol, ycol, leny, iylast, iystop, iystrt;

  /* Parameter adjustments */
  --y;
  --x;
  --xpnt;

  /* Function Body */
  mm = m;
  iylast = 0;
  leny = ldy;
  /*       ------------------------------------ */
  /*       TO COMPUTE EACH COLUMN YCOL OF Y ... */
  /*       ------------------------------------ */
  for (ycol = 1; ycol <= q; ++ycol) {
    iystrt = iylast + 1;
    iystop = iystrt + mm - 1;
    iylast += leny;
    /*           -------------------------------------------------- */
    /*           ... PERFORM THE APPROPRIATE MATRIX VECTOR MULTIPLY: */
    /*               X * A(*,YCOL). */
    /*           -------------------------------------------------- */
    for (xcol = 1; xcol <= n; ++xcol) {
      j1 = xpnt[xcol];
      i1 = xpnt[xcol + 1] - mm;
      z__2 = -x[j1];
      z__1 = z__2 * x[i1];
      a1 = z__1;
      for (iy = iystrt; iy <= iystop; ++iy) {
        i__4 = iy;
        i__5 = iy;
        i__6 = i1;
        z__2 = a1 * x[i__6];
        z__1 = y[i__5] + z__2;
        y[i__4] = z__1;
        ++i1;
      }
    }
    --mm;
    --leny;
  }

  return 0;
} /* mmpy1 */

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

template <typename Scalar, typename Index>
void chlsup(
  Index m, Index n, Index* __restrict split, Index* __restrict xpnt,
  Scalar* __restrict x, Index& iflag) {
  Index q, mm, nn, jblk, jpnt;
  Index fstcol, nxtcol;

  /* Parameter adjustments */
  --x;
  --xpnt;
  --split;

  /* Function Body */
  jblk = 0;
  fstcol = 1;
  mm = m;
  jpnt = xpnt[fstcol];

  /*       ---------------------------------------- */
  /*       FOR EACH BLOCK JBLK IN THE SUPERNODE ... */
  /*       ---------------------------------------- */
  while (fstcol <= n) {
    ++jblk;
    nn = split[jblk];
    /*           ------------------------------------------ */
    /*           ... PERFORM PARTIAL CHOLESKY FACTORIZATION */
    /*               ON THE BLOCK. */
    /*           ------------------------------------------ */
    pchol(mm, nn, &xpnt[fstcol], &x[1], iflag);
    if (iflag == 1) {
      return;
    }
    /*           ---------------------------------------------- */
    /*           ... APPLY THE COLUMNS IN JBLK TO ANY COLUMNS */
    /*               OF THE SUPERNODE REMAINING TO BE COMPUTED. */
    /*           ---------------------------------------------- */
    nxtcol = fstcol + nn;
    q = n - nxtcol + 1;
    mm -= nn;
    jpnt = xpnt[nxtcol];
    if (q > 0) {
      mmpy1(mm, nn, q, &xpnt[fstcol], &x[1], &x[jpnt], mm);
    }
    fstcol = nxtcol;
  }
} /* chlsup */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* **************     MMPY  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*   INPUT PARAMETERS - */
/*       M               -   NUMBER OF ROWS IN X AND IN Y. */
/*       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS */
/*                           IN A. */
/*       Q               -   NUMBER OF COLUMNS IN A AND Y. */
/*       SPLIT(*)        -   BLOCK PARTITIONING OF X. */
/*       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE */
/*                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO */
/*                           USED TO ACCESS THE ROWS OF A. */
/*       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A. */
/*       LDY             -   LENGTH OF FIRST COLUMN OF Y. */

/*   UPDATED PARAMETERS - */
/*       Y(*)            -   ON OUTPUT, Y = Y + AX. */

/* *********************************************************************** */

template <typename Scalar, typename Index>
void mmpy(
  Index m, Index n, Index q, Index* __restrict split, Index* __restrict xpnt,
  Scalar* __restrict x, Scalar* __restrict y, Index ldy) {
  Index nn, blk, fstcol;
  blk = 0;
  fstcol = 0;
  while (fstcol < n) {
    nn = split[blk];
    mmpy1(m, nn, q, &xpnt[fstcol], x, y, ldy);
    fstcol += nn;
    ++blk;
  }
} /* mmpy */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     MMPYI  .... MATRIX-MATRIX MULTIPLY     ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE - */
/*       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA, */
/*       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY */
/*       CODES. */

/*       MATRIX X HAS ONLY 1 COLUMN. */

/*   INPUT PARAMETERS - */
/*       M               -   NUMBER OF ROWS IN X AND IN Y. */
/*       Q               -   NUMBER OF COLUMNS IN A AND Y. */
/*       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE */
/*                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO */
/*                           USED TO ACCESS THE ROWS OF A. */
/*       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A. */
/*       D               -   DIAGONAL ENTRY IN COLUMN X BACK IN CHOLESKY. */
/*       IY(*)           -   IY(COL) POINTS TO THE BEGINNING OF COLUMN */
/*       RELIND(*)       -   RELATIVE INDICES. */

/*   UPDATED PARAMETERS - */
/*       Y(*)            -   ON OUTPUT, Y = Y + AX. */

/* *********************************************************************** */

template <typename Scalar, typename Index>
void mmpyi_(
  Index m, Index q, Index* __restrict xpnt, Scalar* __restrict x, Scalar d,
  Index* __restrict iy, Scalar* __restrict y, Index* __restrict relind) {
  /* Local variables */
  Scalar a, z1, z2;
  Index i, k, col, isub, ylast;

  /* Parameter adjustments */
  --relind;
  --y;
  --iy;
  --x;
  --xpnt;

  /* Function Body */
  for (k = 1; k <= q; ++k) {
    col = xpnt[k];
    ylast = iy[col + 1] - 1;
    z2 = -(d);
    z1 = z2 * x[k];
    a = z1;
    for (i = k; i <= m; ++i) {
      isub = xpnt[i];
      isub = ylast - relind[isub];
      y[isub] += a * x[i];
    }
  }
} /* mmpyi_ */

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

template <typename Scalar, typename Index>
void blkfc2(
  Index nsuper, Index* xsuper, Index* snode, Index* split, Index* xlindx,
  Index* lindx, Index* xlnz, Scalar* lnz, Index* link, Index* length,
  Index* indmap, Index* relind, Index tmpsiz, Scalar* temp, Index& iflag) {
  /* Local variables */
  Index i, ilen, jlen, klen, jsup, ksup;
  Index fjcol, fkcol, ljcol;
  Index klast, kdpnt, ilpnt, jlpnt, klpnt, store;
  Index jxpnt, kxpnt, inddif;
  Index njcols, nkcols;
  Index ncolup, kfirst, nxtcol, nxksup, nxtsup;

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
  iflag = 0;

  /*       ----------------------------------------------------------- */
  /*       INITIALIZE EMPTY ROW LISTS IN LINK(*) AND ZERO OUT TEMP(*). */
  /*       ----------------------------------------------------------- */
  std::fill(link + 1, link + nsuper + 1, Index(0));
  std::fill(temp + 1, temp + tmpsiz + 1, Scalar(0));

  /*       --------------------------- */
  /*       FOR EACH SUPERNODE JSUP ... */
  /*       --------------------------- */
  for (jsup = 1; jsup <= nsuper; ++jsup) {
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
    ldindx(jlen, &lindx[jxpnt], &indmap[1]);

    /*           ----------------------------------------- */
    /*           FOR EVERY SUPERNODE KSUP IN ROW(JSUP) ... */
    /*           ----------------------------------------- */
    ksup = link[jsup];
L300:
    if (ksup > 0) {
      nxksup = link[ksup];

      /*               -------------------------------------------------------
						 */
      /*               GET INFO ABOUT THE CMOD(JSUP,KSUP) UPDATE. */

      /*               FKCOL  ...  FIRST COLUMN OF SUPERNODE KSUP. */
      /*               NKCOLS ...  NUMBER OF COLUMNS IN SUPERNODE KSUP. */
      /*               KLEN   ...  LENGTH OF ACTIVE PORTION OF COLUMN FKCOL. */
      /*               KXPNT  ...  POINTER TO INDEX OF FIRST NONZERO IN ACTIVE
						 */
      /*                           PORTION OF COLUMN FJCOL. */
      /*               -------------------------------------------------------
						 */
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

        for (i = 0; i < klen; ++i) {
          nxtcol = lindx[kxpnt + i];
          if (nxtcol > ljcol) {
            goto L500;
          }
        }
        i = klen;
L500:
        ncolup = i;

        if (nkcols == 1) {
          /*                       ----------------------------------------------
								 */
          /*                       UPDATING TARGET SUPERNODE BY TRIVIAL */
          /*                       SUPERNODE (WITH ONE COLUMN). */

          /*                       KLPNT  ...  POINTER TO FIRST NONZERO IN
								 * ACTIVE */
          /*                                   PORTION OF COLUMN FKCOL. */
          /*                       KDPNT  ...  POINTER TO DIAGONAL ENTRY OF */
          /*                                   COLUMN FKCOL. */
          /*                       ----------------------------------------------
								 */
          klpnt = xlnz[fkcol + 1] - klen;
          kdpnt = xlnz[fkcol];
          mmpyi_(
            klen, ncolup, &lindx[kxpnt], &lnz[klpnt], lnz[kdpnt], &xlnz[1],
            &lnz[1], &indmap[1]);
        } else {
          /*                       --------------------------------------------
								 */
          /*                       KFIRST ...  FIRST INDEX OF ACTIVE PORTION OF
								 */
          /*                                   SUPERNODE KSUP (FIRST COLUMN TO
								 */
          /*                                   BE UPDATED). */
          /*                       KLAST  ...  LAST INDEX OF ACTIVE PORTION OF
								 */
          /*                                   SUPERNODE KSUP. */
          /*                       --------------------------------------------
								 */

          kfirst = lindx[kxpnt];
          klast = lindx[kxpnt + klen - 1];
          inddif = indmap[kfirst] - indmap[klast];

          if (inddif < klen) {
            /*                           ---------------------------------------
									 */
            /*                           DENSE CMOD(JSUP,KSUP). */

            /*                           ILPNT  ...  POINTER TO FIRST NONZERO IN
									 */
            /*                                       COLUMN KFIRST. */
            /*                           ILEN   ...  LENGTH OF COLUMN KFIRST. */
            /*                           ---------------------------------------
									 */
            ilpnt = xlnz[kfirst];
            ilen = xlnz[kfirst + 1] - ilpnt;
            mmpy(
              klen, nkcols, ncolup, &split[fkcol], &xlnz[fkcol], &lnz[1],
              &lnz[ilpnt], ilen);
          } else {
            /*                           ------------------------------- */
            /*                           GENERAL SPARSE CMOD(JSUP,KSUP). */
            /*                           COMPUTE CMOD(JSUP,KSUP) UPDATE */
            /*                           IN WORK STORAGE. */
            /*                           ------------------------------- */
            store = klen * ncolup - ncolup * (ncolup - 1) / 2;
            if (store > tmpsiz) {
              iflag = -2;
              return;
            }
            mmpy(
              klen, nkcols, ncolup, &split[fkcol], &xlnz[fkcol], &lnz[1],
              &temp[1], klen);
            /*                           ----------------------------------------
									 */
            /*                           GATHER INDICES OF KSUP RELATIVE TO
									 * JSUP. */
            /*                           ----------------------------------------
									 */
            igathr(klen, &lindx[kxpnt], &indmap[1], &relind[1]);
            /*                           --------------------------------------
									 */
            /*                           INCORPORATE THE CMOD(JSUP,KSUP) BLOCK
									 */
            /*                           UPDATE INTO THE TO APPROPRIATE COLUMNS
									 */
            /*                           OF L. */
            /*                           --------------------------------------
									 */
            assmb(
              klen, ncolup, &temp[1], &relind[1], &xlnz[fjcol], &lnz[1], jlen);
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
        mmpy(
          klen, nkcols, njcols, &split[fkcol], &xlnz[fkcol], &lnz[1],
          &lnz[jlpnt], jlen);
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
    chlsup(jlen, njcols, &split[fjcol], &xlnz[fjcol], &lnz[1], iflag);
    if (iflag != 0) {
      iflag = -1;
      return;
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
  }
} /* blkfc2 */

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
/*       TMPVEC          -   COMPLEX*16 WORKING STORAGE OF LENGTH */
/*                           NEQNS. */

/* *********************************************************************** */

template <typename Scalar, typename Index>
void blkfct(
  Index neqns, Index nsuper, Index* xsuper, Index* snode, Index* split,
  Index* xlindx, Index* lindx, Index* xlnz, Scalar* lnz, Scalar* diag,
  Index iwsiz, Index* iwork, Index tmpsiz, Scalar* tmpvec, Index& iflag) {
  /* Parameter adjustments */
  --diag;
  --lnz;
  --xlnz;

  /* Function Body */
  iflag = 0;
  if (iwsiz < 2 * neqns + 2 * nsuper) {
    iflag = -3;
    std::cerr << "\n";
    std::cerr << "*** INTEGER WORK SPACE = " << iwsiz << "\n";
    std::cerr << "*** IS SMALLER THAN REQUIRED = " << 2 * nsuper + 2 * neqns
              << "\n";
    std::cerr << "\n";
    return;
  }

  blkfc2(
    nsuper, xsuper, snode, split, xlindx, lindx, &xlnz[1], &lnz[1], iwork,
    &iwork[nsuper], &iwork[2 * nsuper], &iwork[2 * nsuper + neqns], tmpsiz,
    tmpvec, iflag);

  if (iflag == -1) {
    std::cerr << "\n";
    std::cerr << "*** MATRIX IS SINGULAR ***\n";
    std::cerr << "*** ZERO DIAGONAL ENTRY ENCOUNTERED ***\n";
    std::cerr << "\n";
    return;
  } else if (iflag == -2) {
    std::cerr << "\n";
    std::cerr << "*** INSUFFICIENT WORK STORAGE [TMPVEC(*)] ***\n";
    std::cerr << "\n";
    return;
  }
  if (iflag == 0) {
    for (Index jcol = 1; jcol <= neqns; ++jcol) {
      diag[jcol] = lnz[xlnz[jcol]];
    }
  }
} /* blkfct */

}}}  // namespace NgPeytonCpp::details::f2c

#endif  // NGPEYTONCPP_DETAILS_NUMERICFACTORIZATION_H
