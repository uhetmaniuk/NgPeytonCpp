/* assmb.f -- translated by f2c (version 20230428).
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
/* ************     ASSMB .... INDEXED ASSEMBLY OPERATION     ************ */
/* *********************************************************************** */
/* *********************************************************************** */

/*   PURPOSE: */
/*       THIS ROUTINE PERFORMS AN INDEXED ASSEMBLY (I.E., SCATTER-ADD) */
/*       OPERATION, ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE */
/*       CHOLESKY CODES. */

/*   INPUT PARAMETERS: */
/*       M               -   NUMBER OF ROWS IN Y. */
/*       Q               -   NUMBER OF COLUMNS IN Y. */
/*       Y               -   BLOCK UPDATE TO BE INCORPORATED INTO FACTOR */
/*                           STORAGE. */
/*       RELIND          -   RELATIVE INDICES FOR MAPPING THE UPDATES */
/*                           ONTO THE TARGET COLUMNS. */
/*       XLNZ            -   POINTERS TO THE START OF EACH COLUMN IN THE */
/*                           TARGET MATRIX. */

/*   OUTPUT PARAMETERS: */
/*       LNZ             -   CONTAINS COLUMNS MODIFIED BY THE UPDATE */
/*                           MATRIX. */

/* *********************************************************************** */

/* Subroutine */ int assmb_(integer *m, integer *q, doublereal *y, integer *
	relind, integer *xlnz, doublereal *lnz, integer *lda)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer ir, il1, iy1, icol, ycol, lbot1, yoff1;


/* *********************************************************************** */

/*       ----------- */
/*       PARAMETERS. */
/*       ----------- */


/*       ---------------- */
/*       LOCAL VARIABLES. */
/*       ---------------- */


/* *********************************************************************** */


    /* Parameter adjustments */
    --lnz;
    --xlnz;
    --relind;
    --y;

    /* Function Body */
    yoff1 = 0;
    i__1 = *q;
    for (icol = 1; icol <= i__1; ++icol) {
	ycol = *lda - relind[icol];
	lbot1 = xlnz[ycol + 1] - 1;
/* DIR$ IVDEP */
	i__2 = *m;
	for (ir = icol; ir <= i__2; ++ir) {
	    il1 = lbot1 - relind[ir];
	    iy1 = yoff1 + ir;
	    lnz[il1] += y[iy1];
	    y[iy1] = 0.;
/* L100: */
	}
	yoff1 = iy1 - icol;
/* L200: */
    }

    return 0;
} /* assmb_ */

