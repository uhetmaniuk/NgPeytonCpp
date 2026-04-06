#ifndef NGPEYTONCPP_DETAILS_CONVERSION_H
#define NGPEYTONCPP_DETAILS_CONVERSION_H


namespace NgPeytonCpp { namespace details { namespace f2c {

/*********************************************************************/
/*********************************************************************/

/*     FLO2HO */

/*        This routine converts the lower triangular part of symmetric */
/*        sparse matrix to its full representation (i.e. both */
/*        upper and lower triangular parts are stored. The diagonal */
/*        of the matrix is INCLUDED in newvals in this version. */

/*     Arguments: */

/*     N        INTEGER (INPUT) */
/*              Dimension of the matrix. */

/*     COLPTR   INTEGER array of size N+1. (INPUT) */
/*              The column pointers of the lower triangular input matrix. */

/*     ROWIND   INTEGER array of size COLPTR(N+1)-1. (INPUT) */
/*              The row indices of the lower triangular input matrix. */

/*     NZVALS   COMPLEX*16 array of size COLPTR(N+1)-1. (INPUT) */
/*              The nonzero entries in the lower triangular part of the */
/*              matrix column-wise compressed. */

/*     NEWPTR   INTEGER array of size N+1. (OUTPUT) */
/*              The column pointers of the converted matrix. */

/*     NEWIND   INTEGER array of size H_NNZ. (OUTPUT) */
/*              The row indices of the converted matrix. */

/*     NEWVALS  COMPLEX*16 array of size H_NNZ (OUTPUT) */
/*              The nonzero entries of the upper and lower triangular */
/*              part of the matrix column-wise compressed. */
/* The diagonal of the matrix is not included. */

/*     IP       INTEGER array of size N. (WORK) */
/*              Work array for storing pointers. */

/*     === work array to accumulate the non-zero */
/*         counts for each row === */

template <typename Scalar, typename Index>
void flo2ho(
  Index n, const Index* colptr, const Index* rowind, const Scalar* nzvals,
  Index* newptr, Index* newind, Scalar* newvals, Index* ip) {
  /* Local variables */
  Index i, j, iend, inew, irow, ibegin;

  for (i = 0; i < n; ++i)
    ip[i] = 0;

  /* Parameter adjustments */
  --ip;

  --newptr;
  --newind;
  --newvals;

  /*     === nonzero count for each row of the */
  /*         lower triangular part */
  /*         (including the diagonal) === */

  for (j = 0; j < n; ++j) {
    for (i = colptr[j]; i < colptr[j + 1]; ++i) {
      irow = rowind[i - 1];
      ++ip[irow];
    }
  }

  /*     === nonzero count for each column */
  /*         (excluding the diagonal)      === */

  for (j = 0; j < n; ++j) {
    ibegin = colptr[j];
    iend = colptr[j + 1] - 1;
    if (iend > ibegin) {
      ip[j + 1] += iend - ibegin;
    }
  }

  /*     === compute pointers to the beginning of each column === */

  newptr[1] = 1;
  for (i = 1; i <= n; ++i) {
    newptr[i + 1] = newptr[i] + ip[i];
  }

  for (i = 1; i <= n; ++i)
    ip[i] = newptr[i];

  /*     === copy the upper triangular part === */
  /*            (excluding the diagonal) */

  for (j = 0; j < n; ++j) {
    ibegin = colptr[j];
    iend = colptr[j + 1] - 1;
    if (ibegin < iend) {
      for (i = ibegin + 1; i <= iend; ++i) {
        irow = rowind[i - 1];
        newind[ip[irow]] = j + 1;
        newvals[ip[irow]] = nzvals[i - 1];
        ++ip[irow];
      }
    }
  }

  /*     === copy the lower triangular part === */
  /*            (including the diagonal) */
  for (j = 0; j < n; ++j) {
    inew = ip[j + 1];
    for (i = colptr[j]; i < colptr[j + 1]; ++i) {
      newind[inew] = rowind[i - 1];
      newvals[inew] = nzvals[i - 1];
      ++inew;
    }
  }
} /* flo2ho */

/*********************************************************************/
/*********************************************************************/

/*     Purpose: */

/*        This routine converts the lower triangular representation */
/*        of symmetric sparse matrix pattern to its full representation */
/*        (i.e. both upper and lower triangular indices are stored. */
/*        Diagonal entries are EXCLUDED from NEWPTR and NEWIND. */

/*     Arguments: */

/*     N        INTEGER (INPUT) */
/*              Dimension of the matrix. */

/*     HNNZ     INTEGER (OUTPUT) */
/*              The number of nonzeros in the entire (both upper and lower */
/*              triangular part) of the matrix. */

/*     COLPTR   INTEGER array of size N+1. (INPUT) */
/*              The column pointers of the lower triangular input matrix. */

/*     ROWIND   INTEGER array of size COLPTR(N+1)-1. (INPUT) */
/*              The row indices of the lower triangular input matrix. */

/*     NEWPTR   INTEGER array of size N+1. (OUTPUT) */
/*              The column pointers of the converted matrix. */

/*     NEWIND   INTEGER array of size H_NNZ. (OUTPUT) */
/*              The row indices of the converted matrix. */

/*     IP       INTEGER array of size N. (WORK) */
/*              Work array for pointers. */

template <typename Index>
void ilo2ho(
  Index n, Index& hnnz, const Index* colptr, const Index* rowind, Index* newptr,
  Index* newind, Index* ip) {
  /* Local variables */
  Index i, j, nnz, iend, inew, irow, ibegin;

  nnz = colptr[n] - 1;

  /*     === work array to accumulate the non-zero */
  /*         counts for each row === */

  for (i = 0; i < n; ++i) {
    ip[i] = 0;
  }

  /* Parameter adjustments */
  --ip;
  --newptr;
  --colptr;
  --rowind;
  --newind;

  /*     === nonzero count for each row of the */
  /*         lower triangular part === */

  for (j = 1; j <= n; ++j) {
    ibegin = colptr[j];
    iend = colptr[j + 1] - 1;
    for (i = ibegin; i <= iend; ++i) {
      irow = rowind[i];
      ++ip[irow];
    }
  }

  /*     === nonzero count for each column === */

  for (j = 1; j <= n; ++j) {
    ibegin = colptr[j];
    iend = colptr[j + 1] - 1;
    if (iend >= ibegin) {
      ip[j] = ip[j] + iend - ibegin - 1;
    }
  }

  /*     === compute pointers to the beginning of each column === */

  newptr[1] = 1;
  for (i = 1; i <= n; ++i) {
    newptr[i + 1] = newptr[i] + ip[i];
  }
  hnnz = newptr[n + 1] - 1;

  for (i = 1; i <= n; ++i) {
    ip[i] = newptr[i];
  }

  /*     === copy the upper triangular part === */

  for (j = 1; j <= n; ++j) {
    ibegin = colptr[j];
    iend = colptr[j + 1] - 1;
    for (i = ibegin + 1; i <= iend; ++i) {
      irow = rowind[i];
      newind[ip[irow]] = j;
      ++ip[irow];
    }
  }

  /*     === copy the lower triangular part === */

  for (j = 1; j <= n; ++j) {
    ibegin = colptr[j];
    iend = colptr[j + 1] - 1;
    inew = ip[j];
    for (i = ibegin + 1; i <= iend; ++i) {
      newind[inew] = rowind[i];
      ++inew;
    }
  }
} /* ilo2ho */

}}}  // namespace NgPeytonCpp::details::f2c

#endif  // NGPEYTONCPP_DETAILS_CONVERSION_H
