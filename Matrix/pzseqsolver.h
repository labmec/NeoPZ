/* Generated by Together */

#ifndef TPZSEQUENCESOLVER_H
#define TPZSEQUENCESOLVER_H
#include "pzsolve.h"
#include "pzstack.h"

class TPZFMatrix;

/**
   Defines sequence solvers
   @ingroup solver
*/
class TPZSequenceSolver : public TPZMatrixSolver {
public:
  /**
     Constructor with initialization parameter
     @param refmat Sets reference matrix to NILL
  */
  TPZSequenceSolver(TPZMatrix *refmat = 0);
  /**
     Copy constructor
     @param copy Model object to be copied from
  */
  TPZSequenceSolver(const TPZSequenceSolver & copy);
  
  void Solve(const TPZFMatrix &F, TPZFMatrix &result, TPZFMatrix *residual = 0);
  
  void ResetSolver();
  
  void AppendSolver(TPZMatrixSolver & solve);
  
  virtual TPZSolver * Clone() const;

 private:    
  TPZStack < TPZMatrixSolver * > fSolvers;
};
#endif //TPZSEQUENCESOLVER_H
