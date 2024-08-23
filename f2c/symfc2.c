/* symfc2.f -- translated by f2c (version 20230428).
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

/* Subroutine */ int symfc2_(integer *neqns, integer *adjlen, integer *xadj, 
	integer *adjncy, integer *perm, integer *invp, integer *colcnt, 
	integer *nsuper, integer *xsuper, integer *snode, integer *nofsub, 
	integer *xlindx, integer *lindx, integer *xlnz, integer *mrglnk, 
	integer *rchlnk, integer *marker, integer *flag__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, knz, head, node, tail, pcol, newi, jptr, kptr, jsup, 
	    ksup, psup, nzbeg, nzend, width, nexti, point, jnzbeg, knzbeg, 
	    length, jnzend, jwidth, fstcol, knzend, lstcol;


/* *********************************************************************** */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* *********************************************************************** */

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
    *flag__ = 0;
    if (*neqns <= 0) {
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
    tail = *neqns + 1;
    point = 1;
    i__1 = *neqns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	marker[i__] = 0;
	xlnz[i__] = point;
	point += colcnt[i__];
/* L50: */
    }
    xlnz[*neqns + 1] = point;
    point = 1;
    i__1 = *nsuper;
    for (ksup = 1; ksup <= i__1; ++ksup) {
	mrglnk[ksup] = 0;
	fstcol = xsuper[ksup];
	xlindx[ksup] = point;
	point += colcnt[fstcol];
/* L100: */
    }
    xlindx[*nsuper + 1] = point;

/*       --------------------------- */
/*       FOR EACH SUPERNODE KSUP ... */
/*       --------------------------- */
    i__1 = *nsuper;
    for (ksup = 1; ksup <= i__1; ++ksup) {

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
	    i__2 = jnzbeg;
	    for (jptr = jnzend; jptr >= i__2; --jptr) {
		newi = lindx[jptr];
		++knz;
		marker[newi] = ksup;
		rchlnk[newi] = rchlnk[head];
		rchlnk[head] = newi;
/* L200: */
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
		i__2 = jnzend;
		for (jptr = jnzbeg; jptr <= i__2; ++jptr) {
		    newi = lindx[jptr];
L400:
		    i__ = nexti;
		    nexti = rchlnk[i__];
		    if (newi > nexti) {
			goto L400;
		    }
		    if (newi < nexti) {
			++knz;
			rchlnk[i__] = newi;
			rchlnk[newi] = nexti;
			marker[newi] = ksup;
			nexti = newi;
		    }
/* L500: */
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
	    i__2 = knzend;
	    for (kptr = knzbeg; kptr <= i__2; ++kptr) {
		newi = adjncy[kptr];
		newi = invp[newi];
		if (newi > fstcol && marker[newi] != ksup) {
/*                       -------------------------------- */
/*                       POSITION AND INSERT NEWI IN LIST */
/*                       AND MARK IT WITH KCOL. */
/*                       -------------------------------- */
		    nexti = head;
L600:
		    i__ = nexti;
		    nexti = rchlnk[i__];
		    if (newi > nexti) {
			goto L600;
		    }
		    ++knz;
		    rchlnk[i__] = newi;
		    rchlnk[newi] = nexti;
		    marker[newi] = ksup;
		}
/* L700: */
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
	i__ = head;
	i__2 = nzend;
	for (kptr = nzbeg; kptr <= i__2; ++kptr) {
	    i__ = rchlnk[i__];
	    lindx[kptr] = i__;
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

/* L1000: */
    }

    return 0;

/*       ----------------------------------------------- */
/*       INCONSISTENCY IN DATA STRUCTURE WAS DISCOVERED. */
/*       ----------------------------------------------- */
L8000:
    *flag__ = -2;
    return 0;

} /* symfc2_ */

