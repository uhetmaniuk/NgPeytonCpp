/* fnsplt.f -- translated by f2c (version 20230428).
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

/* Subroutine */ int fnsplt_(integer *neqns, integer *nsuper, integer *xsuper,
	 integer *xlindx, integer *cachsz, integer *split)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer kcol, used, ksup, cache, ncols, width, height, curcol, 
	    fstcol, lstcol, nxtblk;


/* *********************************************************************** */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */

/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */

/* ******************************************************************* */

/*       -------------------------------------------- */
/*       COMPUTE THE NUMBER OF 8-BYTE WORDS IN CACHE. */
/*       -------------------------------------------- */
    /* Parameter adjustments */
    --split;
    --xlindx;
    --xsuper;

    /* Function Body */
    if (*cachsz <= 0) {
	cache = 2000000000;
    } else {
	cache = (real) (*cachsz) * 1024.f / 8.f * .9f;
    }

/*       --------------- */
/*       INITIALIZATION. */
/*       --------------- */
    i__1 = *neqns;
    for (kcol = 1; kcol <= i__1; ++kcol) {
	split[kcol] = 0;
/* L100: */
    }

/*       --------------------------- */
/*       FOR EACH SUPERNODE KSUP ... */
/*       --------------------------- */
    i__1 = *nsuper;
    for (ksup = 1; ksup <= i__1; ++ksup) {
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
	    used = height << 1;
	    --height;
	}

/*               -------------------------------------- */
/*               ... WHILE THE CACHE IS NOT FILLED AND */
/*                   THERE ARE COLUMNS OF THE SUPERNODE */
/*                   REMAINING TO BE PROCESSED ... */
/*               -------------------------------------- */
L300:
	if (used + height < cache && curcol < lstcol) {
/*                   -------------------------------- */
/*                   ... ADD ANOTHER COLUMN TO CACHE. */
/*                   -------------------------------- */
	    ++curcol;
	    ++ncols;
	    used += height;
	    --height;
	    goto L300;
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
/* L1000: */
    }

    return 0;
} /* fnsplt_ */

