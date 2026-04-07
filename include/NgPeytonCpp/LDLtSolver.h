#ifndef NGPEYTONCPP_LDLTSOLVER_H
#define NGPEYTONCPP_LDLTSOLVER_H

#include <complex>
#include <memory>
#include <type_traits>
#include <vector>

namespace NgPeytonCpp {

/// Ordering strategy for the symbolic factorization.
enum class Ordering {
  MMD = -1,          ///< Multiple Minimum Degree
  UserProvided = 0,  ///< User-provided permutation (perm_ must be set)
  MetisNodeND = 2,   ///< METIS Node Nested Dissection (requires -DMETIS)
  MetisEdgeND = 3    ///< METIS Edge Nested Dissection (requires -DMETIS)
};

template <typename Scalar, typename Index>
class LDLtSolver {
public:
  static_assert(!std::is_const<Index>::value, "Index must not be const");
  static_assert(!std::is_const<Scalar>::value, "Scalar must not be const");
  static_assert(
    std::is_integral<Index>::value, "Index must be an integral type");
  static_assert(
    std::is_same<Scalar, float>::value || std::is_same<Scalar, double>::value ||
      std::is_same<Scalar, std::complex<float>>::value ||
      std::is_same<Scalar, std::complex<double>>::value,
    "Scalar must be float, double, complex<float>, or complex<double>");

  /// \brief Remove default constructor
  LDLtSolver() = delete;

  /// \brief Constructor
  ///
  /// \param[in] n_      Number of equations
  /// \param[in] colptr_ Column pointer array, length n+1 (0-based)
  /// \param[in] rowind_ Row index array, length colptr[n] (0-based)
  /// \param[in] nzvals_ Nonzero values array
  /// \param[in] order_  Ordering strategy
  /// \param[in] perm_   User-provided permutation array (0-based), or nullptr
  /// \note Input data correspond to a CSC sparse matrix.
  LDLtSolver(
    Index n_, const Index* colptr_, const Index* rowind_, const Scalar* nzvals_,
    Ordering order_, const Index* perm_ = nullptr);

  /// \brief Default destructor
  ~LDLtSolver() = default;

  /// \brief Solve A*x = rhs
  ///
  /// \param[in]  rhs  Right-hand side vector, length n
  /// \param[out] x    Solution vector, length n
  void solve(const Scalar* rhs, Scalar* x);

protected:
  // --- Persistent state (needed by solve) ---
  Index n = 0;
  Index nsuper = 0;
  std::vector<Index> xsuper;
  std::vector<Index> xlindx;
  std::vector<Index> lindx;
  std::vector<Index> xlnz;
  std::unique_ptr<Scalar[]> lnz;
  std::unique_ptr<Scalar[]> diag;
  std::vector<Index> perm;
  std::vector<Index> invp;
  std::unique_ptr<Scalar[]> newrhs;

  // --- Temporaries (freed after factorization) ---
  Index nsub = 0;
  Index nnzl = 0;
  Index tmpsiz = 0;
  bool fullrep = false;
  std::vector<Index> snodes = {};
  std::vector<Index> split = {};
  std::vector<Index> xadj = {};
  std::vector<Index> adj = {};
  std::vector<Index> iwork = {};

  /// Build the full adjacency structure (no diagonal) from CSC input.
  /// Populates xadj, adj, iwork, fullrep. Input arrays are 1-based.
  void buildAdjacency(Index nnz, const Index* cp, const Index* ri, Index& nnza);

  /// Compute fill-reducing ordering and symbolic factorization.
  /// Populates perm, invp, colcnt, nsuper, xsuper, snodes, xlindx, lindx,
  /// xlnz, split, tmpsiz, nnzl, nsub.
  void symbolicAnalysis(Index nnza, Ordering order_, const Index* pm);

  /// Numerical LDL' factorization.
  /// Populates lnz, diag, newrhs. Frees temporaries afterwards.
  void numericalFactorization(
    Index nnz, const Index* cp, const Index* ri, const Scalar* nzvals_);
};

}  // namespace NgPeytonCpp

#endif  // NGPEYTONCPP_LDLTSOLVER_H
