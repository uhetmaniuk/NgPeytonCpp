#ifndef NGPEYTONCPP_LDLTSOLVER_IMPL_H
#define NGPEYTONCPP_LDLTSOLVER_IMPL_H

#include "NgPeytonCpp/LDLtSolver.h"

#include "details/Utilities.h"
#include "details/Conversion.h"
#include "details/EliminationTree.h"
#include "details/Ordering.h"
#include "details/SymbolicFactorization.h"
#include "details/BLAS.h"
#include "details/NumericFactorization.h"
#include "details/Solve.h"

#ifdef METIS
extern void METIS_EdgeND(
  int* n, int* xadj, int* adj, int* numflag, int* options, int* perm,
  int* invp);

extern void METIS_NodeND(
  int* n, int* xadj, int* adj, int* numflag, int* options, int* perm,
  int* invp);
#endif

#endif  // NGPEYTONCPP_LDLTSOLVER_IMPL_H
