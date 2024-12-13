#ifndef TPZHDIVAPPROXCREATOR_H
#define TPZHDIVAPPROXCREATOR_H

#include "TPZApproxCreator.h"

class TPZMaterial;
class TPZMultiphysicsCompMesh;

/// @class TPZHDivApproxCreator
/// @brief Manages the construction of approximation spaces with HDiv spaces
class TPZHDivApproxCreator : public TPZApproxCreator {

protected:
    
    /// Type of HDiv family to be used. Some of the options are HDivStandard, HDivKernel, and HDivConstant
    HDivFamily fHDivFam = HDivFamily::EHDivStandard;

    /// A map relating the duplicated connect and its pair
    std::map<int64_t,int64_t> fConnDuplicated; //Check with Jeferson if this variable is indeed necessary
    
    /// A flag to control if the pressure dofs will be condensed (just for the artificial-compressibility iterative solver)
    bool fCondensePressure = false;
public:
    
    /// Default constructor
    TPZHDivApproxCreator() = default;
    
    /// Constructor with gmesh
    TPZHDivApproxCreator(TPZGeoMesh * gmesh);

    /// Default destructor
    ~TPZHDivApproxCreator();

    /// Get/Set HDiv family
    HDivFamily &HdivFamily(){return fHDivFam;}
    /// Get HDiv family
    const HDivFamily &HdivFamily() const {return fHDivFam;}
    
    /// Driver function. Will create the atomic meshes (HDiv, L2, etc.) and an associate multiphysics mesh
    TPZMultiphysicsCompMesh *CreateApproximationSpace() override;
    
    /// Create interface elements on hybridizes spaces
    /// @param mphys multiphysics compmesh
    void AddInterfaceComputationalElements(TPZMultiphysicsCompMesh *mphys);

    /// Function used for debugging. Prints the elements and their connects
    void PrintMeshElementsConnectInfo(TPZCompMesh* cmesh);
    
    /// Function used for debugging. Prints the meshes in txt format
    /// @param mpcmesh Multiphysics compmesh
    static void PrintAllMeshes(TPZMultiphysicsCompMesh* mpcmesh);
    
    void CreateAtomicMeshes(TPZManVector<TPZCompMesh*,7>& meshvec, int& lagLevelCounter);
    
    void CreateMultiPhysicsMesh(TPZManVector<TPZCompMesh*,7>& meshvec, int& lagLevelCounter, TPZMultiphysicsCompMesh*& cmeshmulti);

    void SetShouldCondensePressure(bool isCondensed){
        if (isCondensed && fProbType != ProblemType::EDarcy) {
            std::cout << "The pressure condensation is only available for Darcy problem at the moment." << std::endl;
            DebugStop();
        }
        fCondensePressure = isCondensed;
    }

    bool ShouldCondensePressure() const {
        return fCondensePressure;
    }

protected:

    /// Checks if the parameters provided do not violate any condition for mesh generation
    void CheckSetupConsistency() override;
        
    /// Groups the elements in data structure to be condensed
    /// @param mcmesh multiphysics compmesh with elements to be condensed
    void GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh) override;
    
    
    /// Fix the side orient of neighboring elements to be +1
    /// @param cmesh cmesh to fix side orient
    void FixSideOrientHydridMesh(TPZCompMesh* cmesh);
    
    /// Creates and HDiv approximation space/cmesh
    TPZCompMesh * CreateHDivSpace();

    /// Creates an L2 approximation space/cmesh
    /// @param pOrder polynomial order of the space
    /// @param lagLevel lagrange multiplier level for factorization
    TPZCompMesh * CreateL2Space(int pOrder,const int lagLevel);
    
    
    /// Creates a constant L2 space (order 0). Used for extra spaces related to fIsRBSpaces flag
    /// @param lagLevel lagrange multiplier level for factorization
    TPZCompMesh * CreateConstantSpace(const int lagLevel);
    
    /// Creates the multiphysics compmesh based on the vector of atomic meshes
    /// @param meshvector vector of atomic meshes. For instance HDiv and L2 for Darcy problem
    TPZMultiphysicsCompMesh * CreateMultiphysicsSpace(const TPZVec<TPZCompMesh*> &meshvector);
    
    /// Creates the rotation space for elasticity problems
    TPZCompMesh * CreateRotationSpace(const int pOrder, const int lagLevel);

        
private:
        
    /// Creates the hybridized mesh for the flux
    /// @param compmesh to be created with hybrid data structure
    void CreateHybridHDivMesh(TPZCompMesh* cmesh);
    
    /// Changes the lag level of the lagrange elements connects so it won't be condensed
    /// @param cmesh L2 compmesh with the lagrange multiplier
    void ChangeLagLevel(TPZCompMesh* cmesh, const int newLagLevel);

    /// @brief If the semi-hybridization is chosen, the following methods allows activate or disable the duplicated connects
    /// @param cmesh 
    void ActivateDuplicatedConnects(TPZCompMesh *cmesh);
        
    /// COMENTAR
    /// @param cmesh COMENTAR
    void PartitionDependMatrix(TPZCompMesh *cmesh);
        
    /// COMENTAR
    /// @param cmesh COMENTAR
    void DisableDuplicatedConnects(TPZCompMesh *cmesh);
    
    /// COMENTAR
    /// @param cmesh COMENTAR
    void GroupDependMatrix(TPZCompMesh *cmesh);
    
    /// Semi hybridizes the compmesh
    /// @param cmesh compmesh to semi hybridize
    void SemiHybridizeDuplConnects(TPZCompMesh *cmesh);
        
    /// Creates a vector with the size of nelements with the group to which each element belongs for condensations
    /// @param cmesh computational mesh
    /// @param elementgroup vector with groups by element
    static void AssociateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup);

    /// Creates a vector with the size of nelements with the group to which each element belongs for condensations
    /// @param cmesh computational mesh
    /// @param elementgroup vector with groups by element
    void AssociateElementsDuplConnects(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup);
    
};


#endif
