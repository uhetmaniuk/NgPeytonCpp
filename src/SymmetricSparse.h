#ifndef CPPSELINV_SYMMETRICSPARSE_H
#define CPPSELINV_SYMMETRICSPARSE_H

#include <memory>
#include <vector>

namespace NgPeytonCpp {

template <typename Scalar>
class SymmetricSparse {
public:
  /// \brief Remove default constructor
  SymmetricSparse() = delete;

  /// \brief Constructor
  ///
  /// \param n_      Number of equations
  /// \param colptr_ Column pointer array, length n+1 (0-based)
  /// \param rowind_ Row index array, length colptr[n] (0-based)
  /// \param order_  Ordering: -1 = MMD, 0 = user-provided (perm_)
  /// \param perm_   User-provided permutation array (0-based), or nullptr
  SymmetricSparse(
    int n_, const int* colptr_, const int* rowind_, int order_,
    const int* perm_ = nullptr);

  /// \brief Default destructor
  ~SymmetricSparse() = default;

  /// \brief Factorize
  ///
  /// \tparam Index
  /// \param colptr_ Column pointer array, length n+1 (0-based)
  /// \param rowind_ Row index array, length colptr[n] (0-based)
  /// \param nzvals_ Nonzero values array
  template <typename Index>
  void ldlTFactorize(
    const Index* colptr_, const Index* rowind_, const Scalar* nzvals_);

  /// \brief Solve
  ///
  /// \param x
  /// \param rhs
  void solve(const Scalar* rhs, Scalar* x);

protected:
  template <typename Value>
  class LDLt_factor {
  public:
    int nsuper = 0;
    int nsub = 0;
    int nnzl = 0;
    std::vector<int> xsuper;
    std::vector<int> snodes;
    std::vector<int> xlindx;
    std::vector<int> lindx;
    std::vector<int> xlnz;
    std::vector<Value> lnz;
    std::vector<Value> diag;
    std::vector<Value> tmat;
    std::vector<int> perm;
    std::vector<int> invp;
    std::vector<int> colcnt;
    int tmpsiz = 0;
    std::vector<int> split;
    std::vector<Scalar> newrhs;
  };

  int n = 0;
  bool fullrep = false;
  std::vector<int> xadj;
  std::vector<int> adj;
  std::vector<Scalar> anz;
  std::vector<int> iwork;

  std::unique_ptr<LDLt_factor<Scalar>> factor = nullptr;

  //
  // Protected member functions
  //
};

}  // namespace NgPeytonCpp

#include "SymmetricSparse_impl.h"

#endif  // CPPSELINV_SYMMETRICSPARSE_H
