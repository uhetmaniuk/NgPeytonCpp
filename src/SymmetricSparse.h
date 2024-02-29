#ifndef CPPSELINV_SYMMETRICSPARSE_H
#define CPPSELINV_SYMMETRICSPARSE_H

#include <memory>
#include <vector>

namespace NgPeytonCpp {

template <typename Scalar> class SymmetricSparse {

public:

  /// \brief Remove default constructor
  SymmetricSparse() = delete;

  /// \brief Constructor
  ///
  /// \param n_
  /// \param colptr_
  /// \param rowind_
  /// \param Lnnz_
  /// \param order_
  /// \param perm_
  SymmetricSparse(int n_, const int *colptr_, const int *rowind_, int &Lnnz_,
                  int order_, const int *perm_);

  /// \brief Default destructor
  ~SymmetricSparse() = default;

  /// \brief Factorize
  ///
  /// \tparam Index
  /// \param colptr_
  /// \param rowind_
  /// \param nzvals_
  template <typename Index>
  void ldlTFactorize(const Index *colptr_, const Index *rowind_,
                     const Scalar *nzvals_);

  /// \brief Solve
  ///
  /// \param x
  /// \param rhs
  void solve(Scalar *x_, const Scalar *rhs_);

protected:

  template <typename Value> class LDLt_factor {
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
  std::vector<Scalar> adiag; // <- Does not appear to be used
  std::vector<int> iwork;

  std::unique_ptr<LDLt_factor<Scalar>> factor = nullptr;

  //
  // Protected member functions
  //

};

}

#include "SymmetricSparse_impl.h"

#endif // CPPSELINV_SYMMETRICSPARSE_H
