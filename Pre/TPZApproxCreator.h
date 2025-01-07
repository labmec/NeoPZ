//
// Created by victor on 21/09/2022.
//

#ifndef PZ_TPZAPPROXCREATOR_H
#define PZ_TPZAPPROXCREATOR_H

#include "TPZEnumApproxFamily.h"
#include "pzstack.h"
#include "pzgmesh.h"

class TPZMaterial;
class TPZMultiphysicsCompMesh;

enum class HybridizationType {ENone, EStandard, EStandardSquared, ESemi};
enum class ProblemType {ENone, EElastic, EDarcy, EStokes};

/// @class TPZApproxCreator
/// @brief Manages the construction of approximation spaces
class TPZApproxCreator {
protected:

    /// Sets what type of hybridization should be used between elements, if any
    HybridizationType fHybridType = HybridizationType::ENone;

    /// Variable controling the problem (Darcy, elasticity, etc.). This may influence the number of compmeshes
    ProblemType fProbType = ProblemType::ENone;

    /// Vector of materials. Will be used to to create the compmeshes
    std::map<int,TPZMaterial *> fMaterialVec;

    /// Pointer to geomesh
    TPZGeoMesh *fGeoMesh = nullptr;

    /// Default porder to be used by all the elements in the cmeshes
    int fDefaultPOrder = 1;

    /// @brief This sets hdiv+/hdiv++ internal/bubble POrder: 0- hdiv, 1- hdiv+, 2- hdiv++
    int fExtraInternalPOrder = 0;

    /// Number of compmeshes which is related to the number of state variables. Default is Darcy that has Flux (Hdiv) and Pressure (L2) cmeshes
    int fNumMeshes = 2;

    /// Boolean telling internal equations should be condensed
    bool fShouldCondense = 1;

    /// Enhances space with constant fields that allow for the condensation of internal degrees of freedom
    bool fIsRBSpaces = false;

    /// Struct with all the data regarding hybridization between elements
    struct HybridizationData{
        /// Matid of element at the border or higher domain element
        int fWrapMatId = -123456; // to check if it has been initialized

        /// Matid of element at the interface between elements
        int fInterfaceMatId;

        /// Matid of lagrange multiplier element
        int fLagrangeMatId;

        /// Material ids generated due second hybridization
        int fSecondInterfaceMatId;

        /// Matid of second lagrange multiplier, useful in double hybridizations
        int fSecondLagrangeMatId;

        /// indicates whether the boundary conditions should be hybridized and how many times it should be
        int fHybridizeBCLevel = 0;
    };

    /// Attribute of struct with all the data regarding hybridization between elements
    HybridizationData fHybridizationData;

    /// Predominant GeoMesh dimension element
    MElementType fElementType = MElementType::ENoType;

public:
    
    /// Default constructor
    TPZApproxCreator() = default;

    /// Constructor with gmesh
    TPZApproxCreator(TPZGeoMesh * gmesh);

    /// Default destructor
    ~TPZApproxCreator() = default;

    /// Inserts a domain material to the vector of materials of this class. These are used to build the cmeshes
    /// @param mat pointer to materials to be added
    int InsertMaterialObject(TPZMaterial * mat);

    /// Inserts a bc material to the vector of materials of this class. These are used to build the cmeshes
    /// @param mat boundary condition material to be added
    int InsertMaterialObject(TPZBndCond * mat);

    /// Get/Set GeoMesh
    TPZGeoMesh *GeoMesh(){return fGeoMesh;}
    /// Get GeoMesh
    const TPZGeoMesh *GeoMesh() const {return fGeoMesh;}

    /// Get/Set Hybridization type
    HybridizationType &HybridType(){return fHybridType;}
    /// Get Hybridization type
    const HybridizationType &HybridType() const {return fHybridType;}

    /// Get/Set Enhanced spaces
    bool &IsRigidBodySpaces(){return fIsRBSpaces;}
    /// Get Enhanced spaces
    const bool &IsRigidBodySpaces() const {return fIsRBSpaces;}
    
    HybridizationData& HybridData() {return fHybridizationData;}
    const HybridizationData& HybridData() const {return fHybridizationData;}

    /// Get/Set Problem type
    ProblemType &ProbType(){return fProbType;}
    /// Get Problem type
    const ProblemType &ProbType() const {return fProbType;}

    /// return number of meshes
    const int NumMeshes() const {return fNumMeshes;}

    /// Set cmeshes default polynomial order
    void SetDefaultOrder(const int ord){
        fDefaultPOrder = ord;
    }
    /// Get cmeshes default polynomial order
    int &GetDefaultOrder(){return fDefaultPOrder;}

    /// Set extra internal porder
    void SetExtraInternalOrder(const int ord){
        fExtraInternalPOrder = ord;
    }

    /// Get extra internal porder
    int &GetExtraInternalOrder(){return fExtraInternalPOrder;}

    /// Set if should condense internal equations
    void SetShouldCondense( bool isCondensed){
        fShouldCondense = isCondensed;
    }

    /// Get if should condense internal equations
    bool GetShouldCondense( bool isCondensed){
        return fShouldCondense;
    }

    /// Driver function. Will create the atomic meshes (HDiv, L2, etc.) and an associate multiphysics mesh. Should be implemented in son classes
    virtual TPZMultiphysicsCompMesh *CreateApproximationSpace() = 0;

    /// Print attributes of class
    void Print(std::ostream &out = std::cout);
    
    /// Sets if should hybridize boundary
    void SetHybridizeBoundary() {
        if(fHybridType == HybridizationType::EStandardSquared) fHybridizationData.fHybridizeBCLevel = 2;
        else if(fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::ESemi) fHybridizationData.fHybridizeBCLevel = 1;
        else{
            std::cout << "ERROR! Please set a hybridization type before calling SetToHybridizeBoundary()" << std::endl;
            DebugStop();
        }
    }

protected:
    
    /// Checks if the hybridization neighborhood is correct
    void CheckNeighborhoodHybridization() const;

    /// Compute Periferal Material ids for additional elements generated during hybridization
    void ComputePeriferalMaterialIds(int base = 10);

    /// Adds the geo els related to the hybridization
    void AddHybridizationGeoElements();

    ///This method checks if the current configuration is valid
    virtual void CheckSetupConsistency() = 0;

    ///Get BC material ids
    std::set<int> GetBCMatIds();

    ///Get material ids from which has the same dimension as the mesh
    std::set<int> GetVolumetricMatIds();

    /// Insert periferal material objects related to geometric. objects created during hybridization (for wrap and lagrange objects)
    void InsertWrapAndLagrangeMaterialObjects(TPZMultiphysicsCompMesh *mphys);

    /// Insert interface periferal material objects related to geometric objects created during hybridization
    void InsertInterfaceMaterialObjects(TPZMultiphysicsCompMesh *mphys);

    /// Groups elements to be condensed
    virtual void GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh) = 0;

    /// Changes the internal order of the volumetric connect of all CompEls in a given CMesh
    void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder) const;

    /// Determine the predominant element type in the GeoMesh
    void SetMeshElementType();
};


#endif //PZ_TPZAPPROXCREATOR_H
