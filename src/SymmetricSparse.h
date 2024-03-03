#ifndef CPPSELINV_SYMMETRICSPARSE_H
#define CPPSELINV_SYMMETRICSPARSE_H

#include <memory>
#include <vector>

namespace NgPeytonCpp {

    template<typename Scalar, typename Index = int>
    class SymmetricSparse {

    public:

        /// \brief Remove default constructor
        SymmetricSparse() = delete;

        /// \brief Constructor
        ///
        /// \param n_
        /// \param colptr_
        /// \param rowind_
        /// \param order_
        /// \param perm_ User-provided permutation array
        SymmetricSparse(Index n_, const Index *colptr_, const Index *rowind_,
                        Index order_, const Index *perm_ = nullptr);

        /// \brief Default destructor
        ~SymmetricSparse() = default;

        /// \brief Factorize
        ///
        /// \tparam Index
        /// \param colptr_
        /// \param rowind_
        /// \param nzvals_
        void ldlTFactorize(const Index *colptr_, const Index *rowind_,
                           const Scalar *nzvals_);

        /// \brief Solve
        ///
        /// \param x
        /// \param rhs
        void solve(const Scalar *rhs, Scalar *x);

    protected:

        template<typename Value, typename Idx>
        class LDLt_factor {
        public:
            Idx nsuper = 0;
            Idx nsub = 0;
            Idx nnzl = 0;
            std::vector<Idx> xsuper;
            std::vector<Idx> snodes;
            std::vector<Idx> xlindx;
            std::vector<Idx> lindx;
            std::vector<Idx> xlnz;
            std::vector<Value> lnz;
            std::vector<Value> diag;
            std::vector<Value> tmat;
            std::vector<Idx> perm;
            std::vector<Idx> invp;
            std::vector<Idx> colcnt;
            Idx tmpsiz = 0;
            std::vector<Idx> split;
            std::vector<Scalar> newrhs;
        };

        Index n = 0;
        bool fullrep = false;
        std::vector<Index> xadj;
        std::vector<Index> adj;
        std::vector<Scalar> anz;
        std::vector<Index> iwork;

        std::unique_ptr<LDLt_factor<Scalar, Index>> factor = nullptr;

        //
        // Protected member functions
        //

    };

}

#include "SymmetricSparse_impl.h"

#endif // CPPSELINV_SYMMETRICSPARSE_H
