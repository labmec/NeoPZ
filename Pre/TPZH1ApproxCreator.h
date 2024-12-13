//
// Created by victor on 21/09/2022.
//

#ifndef PZ_TPZH1APPROXCREATOR_H_H
#define PZ_TPZH1APPROXCREATOR_H_H

#include "TPZEnumApproxFamily.h"
#include "TPZApproxCreator.h"
#include "pzgmesh.h"

/// @class TPZHDivApproxCreator
/// @brief Manages the construction of approximation spaces with HDiv spaces
class TPZH1ApproxCreator : public TPZApproxCreator
{
protected:

    /// Type of H1 family to be used. Currently only standard is available.
    H1Family fH1Fam = H1Family::EH1Standard;

public:

    /// Default constructor
    TPZH1ApproxCreator() = default;

    /// Constructor with gmesh
    TPZH1ApproxCreator(TPZGeoMesh * gmesh);

    /// Default destructor
    ~TPZH1ApproxCreator() = default;

    /// Driver function. Will create the atomic meshes (HDiv, L2, etc.) and an associate multiphysics mesh
    TPZMultiphysicsCompMesh *CreateApproximationSpace() override;

    /// Driver function. Will create the atomic meshes (HDiv, L2, etc.) and an associate multiphysics mesh
    TPZCompMesh *CreateClassicH1ApproximationSpace();

    /// Group and condense computational elements
    void GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh) override;


private:

    /// Create boundary HDiv space/cmesh
    TPZCompMesh *CreateBoundaryHDivSpace();

    /// Creates an L2 approximation space/cmesh
    TPZCompMesh * CreateL2Space();

    /// Creates an constant approximation space/cmesh;
    /// The space has a different connect disposition from L2-spaces, easier to read for some users.
    TPZCompMesh * CreateConstantSpace(const int &lagMult);

    /// Insert material objects into L2-mesh
    void InsertL2MaterialObjects(TPZCompMesh * L2Mesh);

    /// Create interface elements on hybridized spaces
    void AddInterfaceComputationalElements(TPZMultiphysicsCompMesh *mcmesh);

    /// Find neighbouring flux element among equal and lower level computational elements.
    /// @param wrapcel is usually a wrap computational element.
    TPZCompEl *FindHDivNeighbouringElement(TPZCompEl *wrapcel);

    /// Creates the multiphysics compmesh based on the vector of atomic meshes
    /// @param meshvector vector of atomic meshes. For instance HDiv and L2 for Darcy problem
    TPZMultiphysicsCompMesh * CreateMultiphysicsSpace(TPZManVector<TPZCompMesh*> meshvector);

    ///This method checks if the current configuration is valid
    void CheckSetupConsistency() override;

    /// Creates the rotation space for elasticity problems
    TPZCompMesh * CreateRotationSpace(){
        DebugStop();
        return nullptr;
    }

    /// Associate elements with a volumetric element in hybridized spaces;
    /// Groups volumetric, wrap and interface elements in standard hybridizations;
    /// Groups volumetric, wrap, interface, lagrange, second interface and bc elements belonging to HDiv spaces for squared hybridizations
    /// Is useless for non-hybridized meshes
    /// Does not (yet?) support semi hybridizations
    /// @param elementgroup maps some comp. elements to a certain volumetric element to be grouped up;
    void AssociateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup);
};

#endif

