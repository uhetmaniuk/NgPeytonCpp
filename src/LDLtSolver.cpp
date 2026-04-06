#include <complex>

#include "LDLtSolver_impl.h"

namespace NgPeytonCpp {

namespace {

/// Convert a 0-based index array to 1-based (internal convention).
template <typename Index>
std::vector<Index> toOneBased(const Index* src, Index count) {
  std::vector<Index> dst(count);
  for (Index i = 0; i < count; ++i)
    dst[i] = src[i] + 1;
  return dst;
}

}  // namespace

template <typename Scalar, typename Index>
LDLtSolver<Scalar, Index>::LDLtSolver(
  Index n_, const Index* colptr_, const Index* rowind_, Ordering order_,
  const Index* perm_)
  : n(n_) {
  Index nnz = 0, nnza = 0, ibegin = 0, iend = 0, i = 0, iwsiz = 0, iflag = 0;
  Index irow = 0, jcol = 0;

  // Convert 0-based input to 1-based for internal use
  auto colptr1 = toOneBased(colptr_, n + 1);
  auto rowind1 = toOneBased(rowind_, colptr_[n]);
  std::vector<Index> perm1;
  if (perm_)
    perm1 = toOneBased(perm_, n);
  const Index* cp = colptr1.data();
  const Index* ri = rowind1.data();
  const Index* pm = perm_ ? perm1.data() : nullptr;

  nnz = cp[n] - 1;

  // Detect full vs lower-triangular representation
  fullrep = details::isFullRepresentation(n, cp, ri);

  // nnza = number of off-diagonal nonzeros in full representation
  nnza = (fullrep) ? nnz - n : 2 * (nnz - n);
  iwsiz = 7 * n + 3;

  // Allocate adjacency storage (full representation, no diagonal)
  xadj.resize(n + 1);
  adj.resize(nnza + n);
  iwork.resize(iwsiz);

  if (fullrep) {
    // Copy structure, removing diagonal entries
    xadj[0] = 1;
    for (i = 0; i < n; i++) {
      xadj[i + 1] = xadj[i] + (cp[i + 1] - cp[i] - 1);
    }
    if (xadj[n] - 1 != nnza) {
      throw std::runtime_error("Inconsistent matrix structure (nnz mismatch)");
    }

    i = 0;
    for (jcol = 1; jcol <= n; jcol++) {
      ibegin = cp[jcol - 1];
      iend = cp[jcol] - 1;
      for (irow = ibegin; irow <= iend; irow++) {
        if (ri[irow - 1] != jcol) {
          adj[i] = ri[irow - 1];
          i++;
        }
      }
    }
  } else {
    // Convert lower-triangular to full representation
    details::f2c::ilo2ho(n, nnza, cp, ri, &xadj[0], &adj[0], &iwork[0]);
  }

  // Copy adjacency for reordering (ordmmd destroys the input)
  std::vector<Index> adj2(adj), xadj2(xadj);

  // Allocate supernode partition storage
  snodes.resize(n);
  xsuper.resize(n + 1);
  colcnt.resize(n);
  perm.resize(n);
  invp.resize(n);

  if (order_ == Ordering::MMD) {
    // Multiple Minimum Degree
    bool sfiflg = true;  // output from ordmmd, not used
    details::f2c::ordmmd<Index>(
      n, &xadj2[0], &adj2[0], &invp[0], &perm[0], iwsiz, &iwork[0], nnzl, nsub,
      &colcnt[0], nsuper, &xsuper[0], &snodes[0], sfiflg, iflag);
    if (iflag != 0) {
      throw std::runtime_error("ordmmd failed (insufficient workspace)");
    }
  }
#ifdef METIS
  else if (order_ == Ordering::MetisNodeND) {
    // METIS Node Nested Dissection
    Index numflag = 1;
    Index options[8];
    options[0] = 0;
    METIS_NodeND(
      &n, &xadj2[0], &adj2[0], &numflag, options, &perm[0], &invp[0]);
  } else if (order_ == Ordering::MetisEdgeND) {
    // METIS Edge Nested Dissection
    Index numflag = 1;
    Index options[8];
    options[0] = 0;
    METIS_EdgeND(
      &n, &xadj2[0], &adj2[0], &numflag, options, &perm[0], &invp[0]);
  }
#endif
  else if (order_ == Ordering::UserProvided) {
    if (!pm) {
      throw std::runtime_error("UserProvided ordering requires a permutation");
    }
    // Use the input permutation vector (already converted to 1-based)
    for (i = 0; i < n; i++) {
      perm[i] = pm[i];
      invp[perm[i] - 1] = i + 1;
    }
  } else {
    throw std::runtime_error("Unknown ordering strategy");
  }

  // Symbolic factorization (not needed when MMD has been called)
  if (order_ != Ordering::MMD) {
    std::fill(iwork.begin(), iwork.end(), Index(0));
    details::sfinit(
      n, nnza, &xadj[0], &adj[0], &perm[0], &invp[0], &colcnt[0], nnzl, nsub,
      nsuper, &snodes[0], &xsuper[0], iwsiz, &iwork[0], iflag);
    if (iflag != 0) {
      throw std::runtime_error("sfinit failed (insufficient workspace)");
    }
  }

  xlindx.resize(n + 1);
  lindx.resize(nsub);
  xlnz.resize(n + 1);

  std::fill(iwork.begin(), iwork.end(), Index(0));

  details::f2c::symfct(
    n, nnza, &xadj[0], &adj[0], &perm[0], &invp[0], &colcnt[0], nsuper,
    &xsuper[0], &snodes[0], nsub, &xlindx[0], &lindx[0], &xlnz[0], iwsiz,
    &iwork[0], iflag);
  if (iflag != 0) {
    throw std::runtime_error("symfct failed (insufficient workspace)");
  }

  // Prepare for numerical factorization and triangular solve
  split.resize(n);

  Index cachsz_ = 700;
  details::f2c::bfinit(
    n, nsuper, &xsuper[0], &snodes[0], &xlindx[0], &lindx[0], cachsz_, tmpsiz,
    &split[0]);
}

/// \brief Numerical LDL' factorization.
/// \param colptr_ Column pointer array, length n+1 (0-based)
/// \param rowind_ Row index array, length colptr[n] (0-based)
/// \param nzvals_ Nonzero values array
template <typename Scalar, typename Index>
void LDLtSolver<Scalar, Index>::ldlTFactorize(
  const Index* colptr_, const Index* rowind_, const Scalar* nzvals_) {
  Index nnz, iwsiz, iflag;
  Index maxsup, jsup, nnzlplus, supsize;

  Index neqns = n;
  iwsiz = 7 * neqns + 3;

  // Convert 0-based input to 1-based for internal use
  auto colptr1 = toOneBased(colptr_, neqns + 1);
  auto rowind1 = toOneBased(rowind_, colptr_[neqns]);
  const Index* cp = colptr1.data();
  const Index* ri = rowind1.data();

  nnz = cp[neqns] - 1;

  anz.resize(2 * nnz - neqns);

  // Allocate extra space for diagonal extraction (not used in LDL')
  maxsup = 0;
  nnzlplus = nnzl;
  for (jsup = 0; jsup < nsuper; jsup++) {
    supsize = xsuper[jsup + 1] - xsuper[jsup];
    maxsup = std::max(maxsup, supsize);
    nnzlplus = nnzlplus + supsize * (supsize - 1) / 2;
  }

  lnz.resize(nnzlplus);
  tmat.resize(tmpsiz);
  diag.resize(neqns);
  newrhs.resize(neqns);

  if (fullrep) {
    std::copy(cp, cp + neqns + 1, xadj.begin());
    std::copy(ri, ri + nnz, adj.begin());
    std::copy(nzvals_, nzvals_ + nnz, anz.begin());
  } else {
    details::f2c::flo2ho(
      neqns, cp, ri, nzvals_, &xadj[0], &adj[0], &anz[0], &iwork[0]);
  }

  details::f2c::inpnv(
    neqns, &xadj[0], &adj[0], &anz[0], &perm[0], &invp[0], nsuper, &xsuper[0],
    &xlindx[0], &lindx[0], &xlnz[0], &lnz[0], iwsiz, &iwork[0], iflag);

  tmat.assign(tmpsiz, Scalar(0));

  details::f2c::blkfct(
    neqns, nsuper, &xsuper[0], &snodes[0], &split[0], &xlindx[0], &lindx[0],
    &xlnz[0], &lnz[0], &diag[0], iwsiz, &iwork[0], tmpsiz, &tmat[0], iflag);

  // Free temporary storage no longer needed after factorization
  std::vector<Index>().swap(iwork);
  std::vector<Index>().swap(xadj);
  std::vector<Index>().swap(adj);
  std::vector<Scalar>().swap(anz);
  std::vector<Scalar>().swap(tmat);
  std::vector<Index>().swap(split);
  std::vector<Index>().swap(colcnt);
  std::vector<Index>().swap(snodes);
}

template <typename Scalar, typename Index>
void LDLtSolver<Scalar, Index>::solve(const Scalar* rhs, Scalar* x) {
  for (Index i = 0; i < n; i++)
    newrhs[i] = rhs[perm[i] - 1];

  details::f2c::blkslv(
    nsuper, &xsuper[0], &xlindx[0], &lindx[0], &xlnz[0], &lnz[0], &newrhs[0]);

  for (Index i = 0; i < n; i++)
    x[i] = newrhs[invp[i] - 1];
}

}  // namespace NgPeytonCpp

// Explicit template instantiations
template class NgPeytonCpp::LDLtSolver<float, int>;
template class NgPeytonCpp::LDLtSolver<double, int>;
template class NgPeytonCpp::LDLtSolver<std::complex<float>, int>;
template class NgPeytonCpp::LDLtSolver<std::complex<double>, int>;
