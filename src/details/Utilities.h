#ifndef NGPEYTONCPP_DETAILS_UTILITIES_H
#define NGPEYTONCPP_DETAILS_UTILITIES_H

#include <algorithm>
#include <vector>


namespace NgPeytonCpp { namespace details {

/// Detect whether a symmetric sparse matrix stored in 1-based CSC format
/// uses a full representation (both triangular halves present) or only
/// the lower triangular part.
///
/// @param n        Number of columns
/// @param colptr   Column pointer array, length n+1 (1-based, internal)
/// @param rowind   Row index array, length colptr[n]-1 (1-based, internal)
/// @return         true if upper triangular entries are present
template <typename Index>
bool isFullRepresentation(
  Index n, const Index* __restrict colptr, const Index* __restrict rowind) {
  for (Index j = 1; j < n; ++j) {
    // Pointer to the last entry of column j
    auto pjend = colptr[j] - 1;
    if (pjend > colptr[j - 1]) {
      // There is at least one off-diagonal entry; find its row index
      auto irow = rowind[pjend - 1];
      // Get column pointer to the first entry of irow-th column
      auto pibeg = colptr[irow - 1];
      // Upper triangular part is present if first entry is not the diagonal
      return rowind[pibeg - 1] != irow;
    }
  }
  return false;
}

/* *********************************************************************** */
/* *********************************************************************** */

// -----------------------------------------------------------------------
// ASSMB — Indexed assembly (scatter-add) operation.
// Authors: Esmond G. Ng and Barry W. Peyton (Oak Ridge, 1994)
// -----------------------------------------------------------------------
template <typename Scalar, typename Index>
void assmb(
  Index m, Index q, Scalar* __restrict y, Index* __restrict relind,
  Index* __restrict xlnz, Scalar* __restrict lnz, Index lda) {
  Index ir, il1, iy1, icol, ycol, lbot1, yoff1;

  yoff1 = 0;
  for (icol = 1; icol <= q; ++icol) {
    ycol = lda - relind[icol - 1];
    lbot1 = xlnz[ycol] - 1;
    for (ir = icol; ir <= m; ++ir) {
      il1 = lbot1 - relind[ir - 1];
      iy1 = yoff1 + ir;
      lnz[il1 - 1] += y[iy1 - 1];
      y[iy1 - 1] = static_cast<Scalar>(0);
    }
    yoff1 = iy1 - icol;
  }
}

// -----------------------------------------------------------------------
// INVINV — Compose two inverse permutations: invp = invp2(invp),
// and compute the forward permutation perm.
// Author: Joseph W.H. Liu (Oak Ridge, 1985)
// Note: perm may alias invp2.
// -----------------------------------------------------------------------
template <typename Index>
void invinv(Index neqns, Index* invp, const Index* invp2, Index* perm) {
  Index i, node, interm;

  for (i = 0; i < neqns; ++i) {
    interm = invp[i];
    invp[i] = invp2[interm - 1];
  }

  for (i = 0; i < neqns; ++i) {
    node = invp[i] - 1;
    perm[node] = i + 1;
  }
}

// -----------------------------------------------------------------------
// FSUP2 — Build the supernode partition vector xsuper from snode.
// Assumes a postordered elimination tree and output from fsup1.
// Authors: Esmond G. Ng and Barry W. Peyton (Oak Ridge, 1992)
// -----------------------------------------------------------------------
template <typename Index>
void fsup2(Index neqns, Index nsuper, const Index* snode, Index* xsuper) {
  Index kcol, ksup, lstsup;

  // Compute xsuper by scanning columns right to left
  lstsup = nsuper + 1;
  for (kcol = neqns; kcol >= 1; --kcol) {
    ksup = snode[kcol - 1];
    if (ksup != lstsup) {
      xsuper[lstsup - 1] = kcol + 1;
    }
    lstsup = ksup;
  }
  xsuper[0] = 1;
}

// -----------------------------------------------------------------------
// FSUP1 — First pass of supernode detection. Computes nsuper and snode
// from the elimination tree and column counts. Assumes postordering.
// Authors: Esmond G. Ng and Barry W. Peyton (Oak Ridge, 1992)
// -----------------------------------------------------------------------
template <typename Index>
void fsup1(
  Index neqns, const Index* etpar, const Index* colcnt, Index& nofsub,
  Index& nsuper, Index* snode) {
  // Compute the fundamental supernode partition
  nsuper = 1;
  snode[0] = 1;
  nofsub = colcnt[0];
  for (Index kcol = 1; kcol < neqns; ++kcol) {
    if (etpar[kcol - 1] == kcol + 1) {
      if (colcnt[kcol - 1] == colcnt[kcol] + 1) {
        snode[kcol] = nsuper;
        continue;
      }
    }
    ++nsuper;
    snode[kcol] = nsuper;
    nofsub += colcnt[kcol];
  }
}

// -----------------------------------------------------------------------
// IGATHR — Gather relative indices from an index map.
// -----------------------------------------------------------------------
template <typename Index = int64_t>
void igathr(
  Index klen, const Index* lindx, const Index* indmap, Index* relind) {
  --indmap;
  for (Index i = 0; i < klen; ++i) {
    relind[i] = indmap[lindx[i]];
  }
}

// -----------------------------------------------------------------------
// LDINDX — Load index vector: maps global row indices to relative
// positions in the supernode column.
// -----------------------------------------------------------------------
template <typename Index = int64_t>
void ldindx(Index jlen, Index* lindx, Index* indmap) {
  --indmap;
  Index curlen = jlen;
  for (Index j = 0; j < jlen; ++j) {
    Index jsub = lindx[j];
    --curlen;
    indmap[jsub] = curlen;
  }
}

}}  // namespace NgPeytonCpp::details

#endif  // NGPEYTONCPP_DETAILS_UTILITIES_H
