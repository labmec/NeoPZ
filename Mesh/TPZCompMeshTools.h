//
//  TPZCompMeshTools.h
//  PZ
//
//  Created by Philippe Devloo on 6/4/15.
//
//

#ifndef __PZ__TPZCompMeshTools__
#define __PZ__TPZCompMeshTools__

#include <stdio.h>
#include "pzcmesh.h"
#include "pzfunction.h"
#include "TPZRenumbering.h"

/**
 * This namespace concentrates general purpose auxiliary methods for computational meshes
 */
namespace TPZCompMeshTools
{
    
    void AddHDivPyramidRestraints(TPZCompMesh *cmesh);
    
    void ExpandHDivPyramidRestraints(TPZCompMesh *cmesh);
    
    void LoadSolution(TPZCompMesh *cpressure, TPZFunction<STATE> &Forcing);

    /// group the elements joining boundary elements and their neighbours
    void GroupElements(TPZCompMesh *cmesh, std::set<int64_t> elbasis, std::set<int64_t> &grouped);
    
    /// Put the element set into a subcompmesh and make the connects internal
    void PutinSubmeshes(TPZCompMesh *cmesh, std::set<int64_t> &elindices, int64_t &index, bool KeepOneLagrangian);
    
    /// Put the element set into a subcompmesh and make the connects internal
    void PutinSubmeshes(TPZCompMesh *cmesh, std::map<int64_t,std::set<int64_t> >&elindices, std::map<int64_t,int64_t> &indices, bool KeepOneLagrangian);
    
    /// group all embedded elements of the computational mesh
    void GroupElements(TPZCompMesh *cmesh);
    
    /// created condensed elements for the elements that have internal nodes
    void CreatedCondensedElements(TPZCompMesh *cmesh, bool KeepOneLagrangian, bool keepmatrix = true);
    
    /// ungroup all embedded elements of the computational mesh
    void UnGroupElements(TPZCompMesh *cmesh);
    
    /// uncondensed elements for the elements that have internal nodes
    void UnCondensedElements(TPZCompMesh *cmesh);

    /// compute the norm of the difference between two meshes
    /// put the computed error in the element solution
    void ComputeDifferenceNorm(TPZCompMesh *mesh1, TPZCompMesh *mesh2, TPZVec<STATE> &square_errors);

    /// adjust the polynomial orders of the hdiv elements such that the internal order is higher than the sideorders
    void AdjustFluxPolynomialOrders(TPZCompMesh *fluxmesh, int hdivplusplus);

    /// set the pressure order acording to the order of internal connect of the elements of the fluxmesh
    void SetPressureOrders(TPZCompMesh *fluxmesh, TPZCompMesh *pressuremesh);
    
    void OptimizeBandwidth(TPZCompMesh *cmesh);

    /// Print cmesh solution per geometric element
    void PrintSolutionByGeoElement(TPZCompMesh *cmesh, std::ostream &out);
}

#endif /* defined(__PZ__TPZCompMeshTools__) */
