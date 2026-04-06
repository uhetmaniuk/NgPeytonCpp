#ifndef NGPEYTONCPP_DETAILS_SYMBOLICFACTORIZATION_H
#define NGPEYTONCPP_DETAILS_SYMBOLICFACTORIZATION_H

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "details/Utilities.h"
#include "details/EliminationTree.h"

namespace NgPeytonCpp { namespace details {
namespace f2c {

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     FNTSIZ ..... COMPUTE WORK STORAGE SIZE FOR BLKFCT     ****** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES THE SIZE OF THE WORKING STORAGE */
/*       REQUIRED BY BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       SNODE           -   SUPERNODE MEMBERSHIP. */
/*       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE. */

/*   OUTPUT PARAMETERS: */
/*       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT. */

/*       RETURNS SIZE OF TEMP ARRAY USED BY BLKFCT FACTORIZATION ROUTINE. */
/*       NOTE THAT THE VALUE RETURNED IS AN ESTIMATE, THOUGH IT IS USUALLY */
/*       TIGHT. */

/*       ---------------------------------------- */
/*       COMPUTE SIZE OF TEMPORARY STORAGE VECTOR */
/*       NEEDED BY BLKFCT. */
/*       ---------------------------------------- */
/* *********************************************************************** */

template <typename Index>
void fntsiz(
  Index nsuper, const Index* xsuper, const Index* snode, const Index* xlindx,
  const Index* lindx, Index& tmpsiz) {
  /* Local variables */
  Index ii, iend, clen, ksup, bound, ncols, width, tsize, ibegin, length,
    cursup, nxtsup;

  /* Parameter adjustments */
  --lindx;
  --xlindx;
  --snode;
  --xsuper;

  /* Function Body */
  tmpsiz = 0;
  for (ksup = nsuper; ksup >= 1; --ksup) {
    ncols = xsuper[ksup + 1] - xsuper[ksup];
    ibegin = xlindx[ksup] + ncols;
    iend = xlindx[ksup + 1] - 1;
    length = iend - ibegin + 1;
    bound = length * (length + 1) / 2;
    if (bound > tmpsiz) {
      cursup = snode[lindx[ibegin]];
      clen = xlindx[cursup + 1] - xlindx[cursup];
      width = 0;
      for (ii = ibegin; ii <= iend; ++ii) {
        nxtsup = snode[lindx[ii]];
        if (nxtsup == cursup) {
          ++width;
          if (ii == iend) {
            if (clen > length) {
              tsize = length * width - (width - 1) * width / 2;
              tmpsiz = std::max(tsize, tmpsiz);
            }
          }
        } else {
          if (clen > length) {
            tsize = length * width - (width - 1) * width / 2;
            tmpsiz = std::max(tsize, tmpsiz);
          }
          length -= width;
          bound = length * (length + 1) / 2;
          if (bound <= tmpsiz) {
            break;
          }
          width = 1;
          cursup = nxtsup;
          clen = xlindx[cursup + 1] - xlindx[cursup];
        }
      }
    }
  }
} /* fntsiz */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ****     FNSPLT ..... COMPUTE FINE PARTITIONING OF SUPERNODES     ***** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE DETERMINES A FINE PARTITIONING OF SUPERNODES */
/*       WHEN THERE IS A CACHE AVAILABLE ON THE MACHINE.  THE FINE */
/*       PARTITIONING IS CHOSEN SO THAT DATA RE-USE IS MAXIMIZED. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       XLINDX          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           POINTERS IN THE SUPERNODE INDICES. */
/*       CACHSZ          -   CACHE SIZE IN KILO BYTES. */
/*                           IF THERE IS NO CACHE, SET CACHSZ = 0. */

/*   OUTPUT PARAMETERS: */
/*       SPLIT           -   INTEGER ARRAY OF SIZE NEQNS CONTAINING THE */
/*                           FINE PARTITIONING. */

/* *********************************************************************** */

template <typename Index>
void fnsplt(
  Index neqns, Index nsuper, const Index* xsuper, const Index* xlindx,
  Index cachsz, Index* split) {
  /* Local variables */
  Index used, ksup, cache, ncols, width, height, curcol, fstcol, lstcol, nxtblk;

  /*       -------------------------------------------- */
  /*       COMPUTE THE NUMBER OF 8-BYTE WORDS IN CACHE. */
  /*       -------------------------------------------- */
  if (cachsz <= 0) {
    cache = 2000000000;
  } else {
    constexpr double tmpf = 1024.0 / 8 * 9.;
    cache = static_cast<Index>(static_cast<double>(cachsz) * tmpf);
  }

  /*       --------------- */
  /*       INITIALIZATION. */
  /*       --------------- */
  std::fill(split, split + neqns, Index(0));

  /* Parameter adjustments */
  --split;
  --xlindx;
  --xsuper;

  /*       --------------------------- */
  /*       FOR EACH SUPERNODE KSUP ... */
  /*       --------------------------- */
  for (ksup = 1; ksup <= nsuper; ++ksup) {
    /*           ----------------------- */
    /*           ... GET SUPERNODE INFO. */
    /*           ----------------------- */
    height = xlindx[ksup + 1] - xlindx[ksup];
    fstcol = xsuper[ksup];
    lstcol = xsuper[ksup + 1] - 1;
    width = lstcol - fstcol + 1;
    nxtblk = fstcol;
    /*           -------------------------------------- */
    /*           ... UNTIL ALL COLUMNS OF THE SUPERNODE */
    /*               HAVE BEEN PROCESSED ... */
    /*           -------------------------------------- */
    curcol = fstcol - 1;
L200:
    /*               ------------------------------------------- */
    /*               ... PLACE THE FIRST COLUMN(S) IN THE CACHE. */
    /*               ------------------------------------------- */
    ++curcol;
    if (curcol < lstcol) {
      ++curcol;
      ncols = 2;
      used = height * 3 - 1;
      height += -2;
    } else {
      ncols = 1;
      used = height * 2;
      --height;
    }

    /*               -------------------------------------- */
    /*               ... WHILE THE CACHE IS NOT FILLED AND */
    /*                   THERE ARE COLUMNS OF THE SUPERNODE */
    /*                   REMAINING TO BE PROCESSED ... */
    /*               -------------------------------------- */
    while (used + height < cache && curcol < lstcol) {
      /*                   -------------------------------- */
      /*                   ... ADD ANOTHER COLUMN TO CACHE. */
      /*                   -------------------------------- */
      ++curcol;
      ++ncols;
      used += height;
      --height;
    }
    /*               ------------------------------------- */
    /*               ... RECORD THE NUMBER OF COLUMNS THAT */
    /*                   FILLED THE CACHE. */
    /*               ------------------------------------- */
    split[nxtblk] = ncols;
    ++nxtblk;
    /*               -------------------------- */
    /*               ... GO PROCESS NEXT BLOCK. */
    /*               -------------------------- */
    if (curcol < lstcol) {
      goto L200;
    }
  }
} /* fnsplt */

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

template <typename Index>
void fcnthn(
  Index neqns, const Index* xadj, const Index* adjncy, const Index* perm,
  const Index* invp, const Index* etpar, Index* rowcnt, Index* colcnt,
  Index& nlnz, Index* set, Index* prvlf, Index* level, Index* weight,
  Index* fdesc, Index* nchild, Index* prvnbr) {
  /* Local variables */
  Index j, k, lca, temp, xsup, last1, last2, lflag, pleaf, hinbr, jstop, jstrt,
    ifdesc, oldnbr, parent, lownbr;

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

  level[0] = 0;
  for (k = neqns; k >= 1; --k) {
    rowcnt[k] = 1;
    colcnt[k] = 0;
    set[k] = k;
    prvlf[k] = 0;
    level[k] = level[etpar[k]] + 1;
    weight[k] = 1;
    fdesc[k] = k;
    nchild[k] = 0;
    prvnbr[k] = 0;
  }
  nchild[0] = 0;
  fdesc[0] = 0;
  for (k = 1; k <= neqns; ++k) {
    parent = etpar[k];
    weight[parent] = 0;
    ++nchild[parent];
    ifdesc = fdesc[k];
    if (ifdesc < fdesc[parent]) {
      fdesc[parent] = ifdesc;
    }
  }
  /*       ------------------------------------ */
  /*       FOR EACH ``LOW NEIGHBOR'' LOWNBR ... */
  /*       ------------------------------------ */
  for (lownbr = 1; lownbr <= neqns; ++lownbr) {
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
    for (j = jstrt; j <= jstop; ++j) {
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
            /*                  ---------------------------------------- */
            /*     ... ACCUMULATE LOWNBR-->HINBR PATH * LENGTH */
            /*                               IN ROWCNT(HINBR). */
            /*                           -----------------------------------------
									 */
            rowcnt[hinbr] = rowcnt[hinbr] + level[lownbr] - level[hinbr];
          } else {
            /*                           -----------------------------------------
									 */
            /*                           ... OTHERWISE, LCA <-- FIND(PLEAF),
									 * WHICH */
            /*                               IS THE LEAST COMMON ANCESTOR OF
									 * PLEAF */
            /*                               AND LOWNBR. */
            /*                               (PATH HALVING.) */
            /*                           -----------------------------------------
									 */
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
            /*                           -------------------------------------
									 */
            /*                           ACCUMULATE PLEAF-->LCA PATH LENGTH IN
									 */
            /*                           ROWCNT(HINBR). */
            /*                           DECREMENT WEIGHT(LCA). */
            /*                           -------------------------------------
									 */
            rowcnt[hinbr] = rowcnt[hinbr] + level[lownbr] - level[lca];
            --weight[lca];
          }
          /*                       ----------------------------------------------
								 */
          /*                       LOWNBR NOW BECOMES ``PREVIOUS LEAF'' OF
								 * HINBR. */
          /*                       ----------------------------------------------
								 */
          prvlf[hinbr] = lownbr;
          lflag = 1;
        }
        /*                   --------------------------------------------------
							 */
        /*                   LOWNBR NOW BECOMES ``PREVIOUS NEIGHBOR'' OF HINBR.
							 */
        /*                   --------------------------------------------------
							 */
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
  nlnz = 0;
  for (k = 1; k <= neqns; ++k) {
    temp = colcnt[k] + weight[k];
    colcnt[k] = temp;
    nlnz += temp;
    parent = etpar[k];
    if (parent != 0) {
      colcnt[parent] += temp;
    }
  }
} /* fcnthn */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  December 27, 1994 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* ******     BFINIT ..... INITIALIZATION FOR BLOCK FACTORIZATION   ****** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS SUBROUTINE COMPUTES ITEMS NEEDED BY THE LEFT-LOOKING */
/*       BLOCK-TO-BLOCK CHOLESKY FACTORITZATION ROUTINE BLKFCT. */

/*   INPUT PARAMETERS: */
/*       NEQNS           -   NUMBER OF EQUATIONS. */
/*       NSUPER          -   NUMBER OF SUPERNODES. */
/*       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING */
/*                           THE SUPERNODE PARTITIONING. */
/*       SNODE           -   SUPERNODE MEMBERSHIP. */
/*       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE. */
/*       CACHSZ          -   CACHE SIZE (IN KBYTES). */

/*   OUTPUT PARAMETERS: */
/*       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT. */
/*       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT */
/*                           INTO CACHE. */

/* *********************************************************************** */

template <typename Index = int64_t>
void bfinit(
  Index neqns, Index nsuper, Index* xsuper, Index* snode, Index* xlindx,
  Index* lindx, Index cachsz, Index& tmpsiz, Index* split) {
  /*       --------------------------------------------------- */
  /*       DETERMINE FLOATING POINT WORKING SPACE REQUIREMENT. */
  /*       --------------------------------------------------- */
  fntsiz(nsuper, xsuper, snode, xlindx, lindx, tmpsiz);

  /*       ------------------------------- */
  /*       PARTITION SUPERNODES FOR CACHE. */
  /*       ------------------------------- */
  fnsplt(neqns, nsuper, xsuper, xlindx, cachsz, split);
} /* bfinit */

/* *********************************************************************** */
/* *********************************************************************** */

/*   Version:        0.3 */
/*   Last modified:  February 13, 1995 */
/*   Authors:        Esmond G. Ng and Barry W. Peyton */

/*   Mathematical Sciences Section, Oak Ridge National Laboratory */

/* *********************************************************************** */
/* *********************************************************************** */
/* *************     SYMFC2 ..... SYMBOLIC FACTORIZATION    ************** */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS ROUTINE PERFORMS SUPERNODAL SYMBOLIC FACTORIZATION ON A */
/*       REORDERED LINEAR SYSTEM.  IT ASSUMES ACCESS TO THE COLUMNS */
/*       COUNTS, SUPERNODE PARTITION, AND SUPERNODAL ELIMINATION TREE */
/*       ASSOCIATED WITH THE FACTOR MATRIX L. */

/*   INPUT PARAMETERS: */
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

/*   OUTPUT PARAMETERS: */
/*       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS */
/*                           INTO THE SUBSCRIPT VECTOR. */
/*       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE */
/*                           COMPRESSED SUBSCRIPTS. */
/*       (I) XLNZ        -   COLUMN POINTERS FOR L. */
/*       (I) FLAG        -   ERROR FLAG: */
/*                               0 - NO ERROR. */
/*                               1 - INCONSISTANCY IN THE INPUT. */

/*   WORKING PARAMETERS: */
/*       (I) MRGLNK      -   ARRAY OF LENGTH NSUPER, CONTAINING THE */
/*                           CHILDREN OF EACH SUPERNODE AS A LINKED LIST. */
/*       (I) RCHLNK      -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE */
/*                           CURRENT LINKED LIST OF MERGED INDICES (THE */
/*                           "REACH" SET). */
/*       (I) MARKER      -   ARRAY OF LENGTH NEQNS USED TO MARK INDICES */
/*                           AS THEY ARE INTRODUCED INTO EACH SUPERNODE'S */
/*                           INDEX SET. */

/* *********************************************************************** */

template <typename Index = int64_t>
int symfc2(
  Index neqns, Index adjlen, Index* xadj, Index* adjncy, Index* perm,
  Index* invp, Index* colcnt, Index nsuper, Index* xsuper, Index* snode,
  Index nofsub, Index* xlindx, Index* lindx, Index* xlnz, Index* mrglnk,
  Index* rchlnk, Index* marker, Index& flag) {
  (void)adjlen;
  (void)nofsub;
  /* Local variables */
  Index i, knz, head, node, tail, pcol, newi, jptr, kptr, jsup, ksup, psup,
    nzbeg, nzend, width, nexti, point, jnzbeg, knzbeg, length, jnzend, jwidth,
    fstcol, knzend, lstcol;

  /* Parameter adjustments */
  --marker;
  --xlnz;
  --snode;
  --colcnt;
  --invp;
  --perm;
  --xadj;
  --adjncy;
  --mrglnk;
  --xlindx;
  --xsuper;
  --lindx;

  /* Function Body */
  flag = 0;
  if (neqns <= 0) {
    return 0;
  }

  /*       --------------------------------------------------- */
  /*       INITIALIZATIONS ... */
  /*           NZEND  : POINTS TO THE LAST USED SLOT IN LINDX. */
  /*           TAIL   : END OF LIST INDICATOR */
  /*                    (IN RCHLNK(*), NOT MRGLNK(*)). */
  /*           MRGLNK : CREATE EMPTY LISTS. */
  /*           MARKER : "UNMARK" THE INDICES. */
  /*       --------------------------------------------------- */
  nzend = 0;
  head = 0;
  tail = neqns + 1;
  point = 1;
  for (i = 1; i <= neqns; ++i) {
    marker[i] = 0;
    xlnz[i] = point;
    point += colcnt[i];
  }
  xlnz[neqns + 1] = point;
  point = 1;
  for (ksup = 1; ksup <= nsuper; ++ksup) {
    mrglnk[ksup] = 0;
    fstcol = xsuper[ksup];
    xlindx[ksup] = point;
    point += colcnt[fstcol];
  }
  xlindx[nsuper + 1] = point;

  /*       --------------------------- */
  /*       FOR EACH SUPERNODE KSUP ... */
  /*       --------------------------- */
  for (ksup = 1; ksup <= nsuper; ++ksup) {
    /*           --------------------------------------------------------- */
    /*           INITIALIZATIONS ... */
    /*               FSTCOL : FIRST COLUMN OF SUPERNODE KSUP. */
    /*               LSTCOL : LAST COLUMN OF SUPERNODE KSUP. */
    /*               KNZ    : WILL COUNT THE NONZEROS OF L IN COLUMN KCOL. */
    /*               RCHLNK : INITIALIZE EMPTY INDEX LIST FOR KCOL. */
    /*           --------------------------------------------------------- */
    fstcol = xsuper[ksup];
    lstcol = xsuper[ksup + 1] - 1;
    width = lstcol - fstcol + 1;
    length = colcnt[fstcol];
    knz = 0;
    rchlnk[head] = tail;
    jsup = mrglnk[ksup];

    /*           ------------------------------------------------- */
    /*           IF KSUP HAS CHILDREN IN THE SUPERNODAL E-TREE ... */
    /*           ------------------------------------------------- */
    if (jsup > 0) {
      /*               --------------------------------------------- */
      /*               COPY THE INDICES OF THE FIRST CHILD JSUP INTO */
      /*               THE LINKED LIST, AND MARK EACH WITH THE VALUE */
      /*               KSUP. */
      /*               --------------------------------------------- */
      jwidth = xsuper[jsup + 1] - xsuper[jsup];
      jnzbeg = xlindx[jsup] + jwidth;
      jnzend = xlindx[jsup + 1] - 1;
      for (jptr = jnzend; jptr >= jnzbeg; --jptr) {
        newi = lindx[jptr];
        ++knz;
        marker[newi] = ksup;
        rchlnk[newi] = rchlnk[head];
        rchlnk[head] = newi;
      }
      /*               ------------------------------------------ */
      /*               FOR EACH SUBSEQUENT CHILD JSUP OF KSUP ... */
      /*               ------------------------------------------ */
      jsup = mrglnk[jsup];
L300:
      if (jsup != 0 && knz < length) {
        /*                   ---------------------------------------- */
        /*                   MERGE THE INDICES OF JSUP INTO THE LIST, */
        /*                   AND MARK NEW INDICES WITH VALUE KSUP. */
        /*                   ---------------------------------------- */
        jwidth = xsuper[jsup + 1] - xsuper[jsup];
        jnzbeg = xlindx[jsup] + jwidth;
        jnzend = xlindx[jsup + 1] - 1;
        nexti = head;
        for (jptr = jnzbeg; jptr <= jnzend; ++jptr) {
          newi = lindx[jptr];
L400:
          i = nexti;
          nexti = rchlnk[i];
          if (newi > nexti) {
            goto L400;
          }
          if (newi < nexti) {
            ++knz;
            rchlnk[i] = newi;
            rchlnk[newi] = nexti;
            marker[newi] = ksup;
            nexti = newi;
          }
        }
        jsup = mrglnk[jsup];
        goto L300;
      }
    }
    /*           --------------------------------------------------- */
    /*           STRUCTURE OF A(*,FSTCOL) HAS NOT BEEN EXAMINED YET. */
    /*           "SORT" ITS STRUCTURE INTO THE LINKED LIST, */
    /*           INSERTING ONLY THOSE INDICES NOT ALREADY IN THE */
    /*           LIST. */
    /*           --------------------------------------------------- */
    if (knz < length) {
      node = perm[fstcol];
      knzbeg = xadj[node];
      knzend = xadj[node + 1] - 1;
      for (kptr = knzbeg; kptr <= knzend; ++kptr) {
        newi = adjncy[kptr];
        newi = invp[newi];
        if (newi > fstcol && marker[newi] != ksup) {
          /*                       -------------------------------- */
          /*                       POSITION AND INSERT NEWI IN LIST */
          /*                       AND MARK IT WITH KCOL. */
          /*                       -------------------------------- */
          nexti = head;
L600:
          i = nexti;
          nexti = rchlnk[i];
          if (newi > nexti) {
            goto L600;
          }
          ++knz;
          rchlnk[i] = newi;
          rchlnk[newi] = nexti;
          marker[newi] = ksup;
        }
      }
    }
    /*           ------------------------------------------------------------ */
    /*           IF KSUP HAS NO CHILDREN, INSERT FSTCOL INTO THE LINKED LIST. */
    /*           ------------------------------------------------------------ */
    if (rchlnk[head] != fstcol) {
      rchlnk[fstcol] = rchlnk[head];
      rchlnk[head] = fstcol;
      ++knz;
    }

    /*           -------------------------------------------- */
    /*           COPY INDICES FROM LINKED LIST INTO LINDX(*). */
    /*           -------------------------------------------- */
    nzbeg = nzend + 1;
    nzend += knz;
    if (nzend + 1 != xlindx[ksup + 1]) {
      goto L8000;
    }
    i = head;
    for (kptr = nzbeg; kptr <= nzend; ++kptr) {
      i = rchlnk[i];
      lindx[kptr] = i;
      /* L800: */
    }

    /*           --------------------------------------------------- */
    /*           IF KSUP HAS A PARENT, INSERT KSUP INTO ITS PARENT'S */
    /*           "MERGE" LIST. */
    /*           --------------------------------------------------- */
    if (length > width) {
      pcol = lindx[xlindx[ksup] + width];
      psup = snode[pcol];
      mrglnk[ksup] = mrglnk[psup];
      mrglnk[psup] = ksup;
    }
  }

  return 0;

  /*       ----------------------------------------------- */
  /*       INCONSISTENCY IN DATA STRUCTURE WAS DISCOVERED. */
  /*       ----------------------------------------------- */
L8000:
  flag = -2;
  return 0;
} /* symfc2 */

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

template <typename Index = int64_t>
void symfct(
  Index neqns, Index adjlen, Index* xadj, Index* adjncy, Index* perm,
  Index* invp, Index* colcnt, Index nsuper, Index* xsuper, Index* snode,
  Index nofsub, Index* xlindx, Index* lindx, Index* xlnz, Index iwsiz,
  Index* iwork, Index& flag) {
  flag = 0;
  if (iwsiz < nsuper + 2 * neqns + 1) {
    flag = -1;
    std::cerr << "\n";
    std::cerr << "*** INTEGER WORK STORAGE = " << iwsiz << "\n";
    std::cerr << "*** IS SMALLER THAN REQUIRED = " << nsuper + 2 * neqns + 1
              << "\n";
    std::cerr << "\n";
    return;
  }
  symfc2(
    neqns, adjlen, xadj, adjncy, perm, invp, colcnt, nsuper, xsuper, snode,
    nofsub, xlindx, lindx, xlnz, iwork, &iwork[nsuper],
    &iwork[nsuper + neqns + 1], flag);
  if (flag == -2) {
    std::cerr << "\n";
    std::cerr << "*** INCONSISTENCY IN THE INPUT\n";
    std::cerr << "*** TO SYMBOLIC FACTORIZATION\n";
    std::cerr << "\n";
    return;
  }
} /* symfct */

}  // namespace f2c

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

template <typename Index>
void sfinit(
  Index neqns, Index nnza, Index* xadj, Index* adjncy, Index* perm, Index* invp,
  Index* colcnt, Index& nnzl, Index& nsub, Index& nsuper, Index* snode,
  Index* xsuper, Index iwsiz, Index* iwork, Index& iflag) {
  (void)nnza;
  iflag = 0;
  if (iwsiz < neqns * 7 + 3) {
    iflag = -1;
    std::cerr << "\n";
    std::cerr << "*** INTEGER WORK STORAGE = " << iwsiz << "\n";
    std::cerr << "*** IS SMALLER THAN REQUIRED = " << 3 + 7 * (neqns) << "\n";
    std::cerr << "\n";
    return;
  }

  /*       ------------------------------------------ */
  /*       COMPUTE ELIMINATION TREE AND POSTORDERING. */
  /*       ------------------------------------------ */
  f2c::etordr(
    neqns, xadj, adjncy, perm, invp, iwork, &iwork[neqns], &iwork[2 * neqns],
    &iwork[3 * neqns]);

  /*       --------------------------------------------- */
  /*       COMPUTE ROW AND COLUMN FACTOR NONZERO COUNTS. */
  /*       --------------------------------------------- */
  f2c::fcnthn(
    neqns, xadj, adjncy, perm, invp, iwork, snode, colcnt, nnzl, &iwork[neqns],
    &iwork[2 * neqns], xsuper, &iwork[3 * neqns], &iwork[4 * neqns + 1],
    &iwork[5 * neqns + 2], &iwork[6 * neqns + 3]);

  /*       --------------------------------------------------------- */
  /*       REARRANGE CHILDREN SO THAT THE LAST CHILD HAS THE MAXIMUM */
  /*       NUMBER OF NONZEROS IN ITS COLUMN OF L. */
  /*       --------------------------------------------------------- */
  f2c::chordr(
    neqns, xadj, adjncy, perm, invp, colcnt, iwork, &iwork[neqns],
    &iwork[2 * neqns], &iwork[3 * neqns]);

  /*       ---------------- */
  /*       FIND SUPERNODES. */
  /*       ---------------- */
  fsup1(neqns, iwork, colcnt, nsub, nsuper, snode);
  fsup2(neqns, nsuper, snode, xsuper);
} /* sfinit */

}  // namespace details
}  // namespace NgPeytonCpp

#endif  // NGPEYTONCPP_DETAILS_SYMBOLICFACTORIZATION_H
