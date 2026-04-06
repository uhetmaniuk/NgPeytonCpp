#ifndef NGPEYTONCPP_LDLTSOLVER_H
#define NGPEYTONCPP_LDLTSOLVER_H

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

  /// \brief Remove default constructor
  LDLtSolver() = delete;

  /// \brief Constructor
  ///
  /// \param n_      Number of equations
  /// \param colptr_ Column pointer array, length n+1 (0-based)
  /// \param rowind_ Row index array, length colptr[n] (0-based)
  /// \param order_  Ordering strategy
  /// \param perm_   User-provided permutation array (0-based), or nullptr
  LDLtSolver(
    Index n_, const Index* colptr_, const Index* rowind_, Ordering order_,
    const Index* perm_ = nullptr);

  /// \brief Default destructor
  ~LDLtSolver() = default;

  /// \brief Factorize
  ///
  /// \param colptr_ Column pointer array, length n+1 (0-based)
  /// \param rowind_ Row index array, length colptr[n] (0-based)
  /// \param nzvals_ Nonzero values array
  void ldlTFactorize(
    const Index* colptr_, const Index* rowind_, const Scalar* nzvals_);

  /// \brief Solve
  ///
  /// \param x
  /// \param rhs
  void solve(const Scalar* rhs, Scalar* x);

protected:
  Index nsuper = 0;
  Index nsub = 0;
  Index nnzl = 0;
  std::vector<Index> xsuper;
  std::vector<Index> snodes;
  std::vector<Index> xlindx;
  std::vector<Index> lindx;
  std::vector<Index> xlnz;
  std::vector<Scalar> lnz;
  std::vector<Scalar> diag;
  std::vector<Scalar> tmat;
  std::vector<Index> perm;
  std::vector<Index> invp;
  std::vector<Index> colcnt;
  Index tmpsiz = 0;
  std::vector<Index> split;
  std::vector<Scalar> newrhs;

  Index n = 0;
  bool fullrep = false;
  std::vector<Index> xadj = {};
  std::vector<Index> adj = {};
  std::vector<Scalar> anz = {};
  std::vector<Index> iwork = {};
};

}  // namespace NgPeytonCpp

#endif  // NGPEYTONCPP_LDLTSOLVER_H
