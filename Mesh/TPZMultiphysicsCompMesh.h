//
//  TPZMultiphysicsCompMesh.h
//  pz
//
//  Created by Omar Durán on 3/21/19.
//

#ifndef TPZMultiphysicsCompMesh_h
#define TPZMultiphysicsCompMesh_h

#include <stdio.h>
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "TPZMultiphysicsInterfaceEl.h"

class TPZMultiphysicsCompMesh : public TPZCompMesh {
    
    /// Vector of active physics: index vector
    /// Define wich space will be active in order to generate equations. Should be defined for each space that you want to use by 0: no active or 1: active
    ///The size have to be the same as the m_mesh_vector
    TPZManVector<MSpaceConfig, 7> m_active_approx_spaces;
    
    /// Vector of computational meshes
    TPZManVector<TPZCompMesh * , 7> m_mesh_vector;
    
public:
    
    /// Default constructor
    TPZMultiphysicsCompMesh();
    
    /// Constructor based on TPZGeoMesh pointer and vector of meshes
    TPZMultiphysicsCompMesh(TPZGeoMesh * gmesh, bool isComplex=false);
    
    /// Constructor based on TPZGeoMesh pointer and vector of meshes
    TPZMultiphysicsCompMesh(TPZAutoPointer<TPZGeoMesh>  gmesh, bool isComplex=false);
    
    /// Copy constructor
    TPZMultiphysicsCompMesh(const TPZMultiphysicsCompMesh &other);
    
    /// Assignement constructor
    TPZMultiphysicsCompMesh & operator=(const TPZMultiphysicsCompMesh &other);
    
    /// Destructor
    ~TPZMultiphysicsCompMesh();
    
    /// Automatic builder for the computational mesh structure
    void AutoBuild() override;
    
    /// Set active approximation spaces
    // active_approx_spaces : vector of the size of mesh_vector containing value 0 or 1
    void BuildMultiphysicsSpace(TPZVec<MSpaceConfig> & active_approx_spaces, const TPZVec<TPZCompMesh * > & mesh_vector);
    
    /// Set active approximation spaces
    void BuildMultiphysicsSpace(TPZVec<TPZCompMesh * > & mesh_vector);
    
    /// Set active approximation spaces
    void BuildMultiphysicsSpace(TPZVec<TPZCompMesh * > & mesh_vector, const TPZVec<int64_t> &gelindexes);
    
    /// Set active approximation spaces
    void BuildMultiphysicsSpaceWithMemory(TPZVec<MSpaceConfig> & active_approx_spaces, TPZVec<TPZCompMesh * > & mesh_vector);
    void BuildMultiphysicsSpaceWithMemory(TPZVec<MSpaceConfig> & active_approx_spaces, TPZVec<TPZCompMesh * > & mesh_vector, std::set<int> matsIdWithMem, std::set<int> matsIdNoMem);
    
    void LoadSolutionFromMeshes();
    
    void LoadSolutionFromMultiPhysics();
    
    /// Get the vector of computational meshes
    TPZVec<TPZCompMesh *> & MeshVector();
    
    /// Get the vector of active physics
    TPZVec<MSpaceConfig> & GetActiveApproximationSpaces();
    
private:
    /// add the elements from the atomic meshes to the multiphysics elements
    void AddElements();
    /// add the connects from the atomic meshes
    void AddConnects();
    
    /// delete the elements and connects
    void CleanElementsConnects();

    template<class TVar>
    void LoadSolutionFromMeshesInternal();
    template<class TVar>
    void LoadSolutionFromMultiPhysicsInternal();
};

#endif /* TPZMultiphysicsCompMesh_h */
