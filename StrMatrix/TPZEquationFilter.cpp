
#include "TPZEquationFilter.h"

void TPZEquationFilter::Scatter(const TPZBaseMatrix &vsmall,
                                TPZBaseMatrix &vexpand) const {
  auto *vsmall_cast_state = dynamic_cast<const TPZFMatrix<STATE> *>(&vsmall);
  auto *vexpand_cast_state = dynamic_cast<TPZFMatrix<STATE> *>(&vexpand);
  if (vsmall_cast_state && vexpand_cast_state) {
    ScatterInternal(*vsmall_cast_state, *vexpand_cast_state);
  }
  // else{
  //   auto *vsmall_cast_cstate =
  //       dynamic_cast<const TPZFMatrix<CSTATE> *>(&vsmall);
  //   auto *vexpand_cast_cstate =
  //       dynamic_cast<TPZFMatrix<CSTATE> *>(&vexpand);
  //   if (vsmall_cast_cstate && vexpand_cast_cstate) {
  //     ScatterInternal(*vsmall_cast_cstate, *vexpand_cast_cstate);
  //   }
  // }
  else {
    PZError << __PRETTY_FUNCTION__;
    PZError << " Incompatible types. Aborting...\n";
    DebugStop();
  }
}

void TPZEquationFilter::Gather(const TPZBaseMatrix &large, TPZBaseMatrix &gathered) const {
  auto *large_cast_state = dynamic_cast<const TPZFMatrix<STATE> *>(&large);
  auto *gathered_cast_state = dynamic_cast<TPZFMatrix<STATE> *>(&gathered);
  if (large_cast_state && gathered_cast_state) {
    GatherInternal(*large_cast_state, *gathered_cast_state);
  }
  // else{
  //   auto *large_cast_cstate =
  //       dynamic_cast<const TPZFMatrix<CSTATE> *>(&large);
  //   auto *gathered_cast_cstate =
  //       dynamic_cast<TPZFMatrix<CSTATE> *>(&gathered);
  //   if (large_cast_cstate && gathered_cast_cstate) {
  //     GatherInternal(*large_cast_cstate, *gathered_cast_cstate);
  //   }
  // }
  else {
    PZError << __PRETTY_FUNCTION__;
    PZError << " Incompatible types. Aborting...\n";
    DebugStop();
  }
}

template <class T>
void TPZEquationFilter::GatherInternal(const TPZFMatrix<T> &large,
                                       TPZFMatrix<T> &gathered) const {
  int64_t neqcondense = this->NActiveEquations();
  if (gathered.Rows() != neqcondense || large.Rows() != fNumEq) {
    DebugStop();
  }
  if (!IsActive()) {
    gathered = large;
    return;
  }
  gathered.Zero();
  for (int64_t i = 0; i < neqcondense; i++)
    gathered(i, 0) = large.GetVal(fActiveEqs[i], 0);
}

template <class TVar>
void TPZEquationFilter::ScatterInternal(const TPZFMatrix<TVar> &vsmall,
                     TPZFMatrix<TVar> &vexpand) const {
  int64_t neqcondense = this->NActiveEquations();
  if (vsmall.Rows() != neqcondense || vexpand.Rows() != fNumEq) {
    DebugStop();
  }
  if (!IsActive()) {
    vexpand = vsmall;
    return;
  }
  vexpand.Zero();

#ifdef PZDEBUG
  {
    for (int64_t i = 0; i < neqcondense; i++) {
      if (fActiveEqs[i] >= fNumEq) {
        DebugStop();
      }
    }
  }
#endif
  for (int64_t i = 0; i < neqcondense; i++)
    vexpand(fActiveEqs[i], 0) = vsmall.GetVal(i, 0);
}

template class TPZRestoreClass<TPZEquationFilter>;
