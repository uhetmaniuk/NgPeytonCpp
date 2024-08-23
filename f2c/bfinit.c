/* bfinit.f -- translated by f2c (version 20230428).
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

/* Subroutine */ int bfinit_(integer *neqns, integer *nsuper, integer *xsuper,
	 integer *snode, integer *xlindx, integer *lindx, integer *cachsz, 
	integer *tmpsiz, integer *split)
{
    extern /* Subroutine */ int fnsplt_(integer *, integer *, integer *, 
	    integer *, integer *, integer *), fntsiz_(integer *, integer *, 
	    integer *, integer *, integer *, integer *);


/* *********************************************************************** */


/* *********************************************************************** */

/*       --------------------------------------------------- */
/*       DETERMINE FLOATING POINT WORKING SPACE REQUIREMENT. */
/*       --------------------------------------------------- */
    /* Parameter adjustments */
    --split;
    --lindx;
    --xlindx;
    --snode;
    --xsuper;

    /* Function Body */
    fntsiz_(nsuper, &xsuper[1], &snode[1], &xlindx[1], &lindx[1], tmpsiz);

/*       ------------------------------- */
/*       PARTITION SUPERNODES FOR CACHE. */
/*       ------------------------------- */
    fnsplt_(neqns, nsuper, &xsuper[1], &xlindx[1], cachsz, &split[1]);

    return 0;
} /* bfinit_ */

