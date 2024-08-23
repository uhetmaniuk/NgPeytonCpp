/* pchol.f -- translated by f2c (version 20230428).
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

/* Subroutine */ int pchol_(integer *m, integer *n, integer *xpnt, doublereal 
	*x, integer *iflag, int (*smxpy)(integer *, integer *, doublereal *, integer *,
            doublereal *))
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer mm;
    static doublereal diag;
    static integer jcol, jpnt;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *);


/* *********************************************************************** */

/*     ----------- */
/*     PARAMETERS. */
/*     ----------- */





/*     ---------------- */
/*     LOCAL VARIABLES. */
/*     ---------------- */



/* *********************************************************************** */

/*       ------------------------------------------ */
/*       FOR EVERY COLUMN JCOL IN THE SUPERNODE ... */
/*       ------------------------------------------ */
    /* Parameter adjustments */
    --x;
    --xpnt;

    /* Function Body */
    mm = *m;
    jpnt = xpnt[1];
    i__1 = *n;
    for (jcol = 1; jcol <= i__1; ++jcol) {

/*           ---------------------------------- */
/*           UPDATE JCOL WITH PREVIOUS COLUMNS. */
/*           ---------------------------------- */
	if (jcol > 1) {
	    i__2 = jcol - 1;
	    (*smxpy)(&mm, &i__2, &x[jpnt], &xpnt[1], &x[1]);
	}

/*           --------------------------- */
/*           COMPUTE THE DIAGONAL ENTRY. */
/*           --------------------------- */
	diag = x[jpnt];
	if (diag == 0.) {
	    *iflag = 1;
	    return 0;
	}
	diag = 1. / diag;

/*           ---------------------------------------------------- */
/*           SCALE COLUMN JCOL WITH RECIPROCAL OF DIAGONAL ENTRY. */
/*           ---------------------------------------------------- */
	--mm;
	++jpnt;
	dscal_(&mm, &diag, &x[jpnt]);
	jpnt += mm;

/* L100: */
    }

    return 0;
} /* pchol_ */

