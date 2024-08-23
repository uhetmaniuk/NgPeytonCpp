/* flo2ho.f -- translated by f2c (version 20230428).
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

/* Subroutine */ int flo2ho_(integer *n, integer *colptr, integer *rowind, 
	doublereal *nzvals, integer *newptr, integer *newind, doublereal *
	newvals, integer *ip)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, nnz, iend, inew, irow, ibegin;


/*     Purpose: */

/*        This routine converts the lower triangluar part of symmetric */
/*        sparse matrix to its full representation (i.e. both */
/*        upper and lower triangular parts are stored.  The diagonal */
/*        of the matrix is INCLUDED in newvals in this version. */

/*     Arguments: */

/*     N        INTEGER (INPUT) */
/*              Dimension of the matrix. */

/*     COLPTR   INTEGER array of size N+1. (INPUT) */
/*              The column pointers of the lower triangular input matrix. */

/*     ROWIND   INTEGER array of size COLPTR(N+1)-1. (INPUT) */
/*              The row indices of the lower triangular input matrix. */

/*     NZVALS   DOUBLE PRECISION array of size COLPTR(N+1)-1. (INPUT) */
/*              The nonzero entries in the lower triangular part of the */
/*              matrix columnwise compressed. */

/*     NEWPTR   INTEGER array of size N+1. (OUTPUT) */
/*              The column pointers of the converted matrix. */

/*     NEWIND   INTEGER array of size H_NNZ. (OUTPUT) */
/*              The row indices of the converted matrix. */

/*     NEWVALS  DOUBLE PRECISION array of size H_NNZ (OUTPUT) */
/*              The nonzero entries of the upper and lower triangular */
/*              part of the matrix columnwize compressed.  The diagonal */
/*              of the matrix is not included. */

/*     IP       INTEGER array of size N. (WORK) */
/*              Work array for storing pointers. */



    /* Parameter adjustments */
    --ip;
    --newptr;
    --colptr;
    --rowind;
    --nzvals;
    --newind;
    --newvals;

    /* Function Body */
    nnz = newptr[*n + 1] - 1;

/*     === work array to accumulate the non-zero */
/*         counts for each row === */

/* vdir noaltcode */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ip[i__] = 0;
    }

/*     === nonzero count for each row of the */
/*         lower triangular part */
/*         (including the diagonal) === */


    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ibegin = colptr[j];
	iend = colptr[j + 1] - 1;
/* vdir nodep */
	i__2 = iend;
	for (i__ = ibegin; i__ <= i__2; ++i__) {
	    irow = rowind[i__];
	    ++ip[irow];
	}
    }

/*     === nonzero count for each column */
/*         (excluding the diagonal)      === */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ibegin = colptr[j];
	iend = colptr[j + 1] - 1;
	if (iend > ibegin) {
	    ip[j] = ip[j] + iend - ibegin;
	}
    }

/*     === compute pointers to the beginning of each column === */

    newptr[1] = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	newptr[i__ + 1] = newptr[i__] + ip[i__];
    }

/* vdir noaltcode */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ip[i__] = newptr[i__];
    }

/*     === copy the upper triangular part === */
/*            (excluding the diagonal) */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ibegin = colptr[j];
	iend = colptr[j + 1] - 1;
	if (ibegin < iend) {
/* vdir nodep */
	    i__2 = iend;
	    for (i__ = ibegin + 1; i__ <= i__2; ++i__) {
		irow = rowind[i__];
		newind[ip[irow]] = j;
		newvals[ip[irow]] = nzvals[i__];
		++ip[irow];
	    }
	}
    }

/*     === copy the lower triangular part === */
/*            (including the diagonal) */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ibegin = colptr[j];
	iend = colptr[j + 1] - 1;
	inew = ip[j];
	if (ibegin <= iend) {
	    i__2 = iend;
	    for (i__ = ibegin; i__ <= i__2; ++i__) {
		newind[inew] = rowind[i__];
		newvals[inew] = nzvals[i__];
		++inew;
	    }
	}
    }
    return 0;
} /* flo2ho_ */

