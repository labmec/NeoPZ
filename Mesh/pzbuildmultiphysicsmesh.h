/*
 *  pzbuildmultiphysicsmesh.h
 *  PZ
 *
 *  Created by Agnaldo on 10/31/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PZBUILDMULTIPHYSICSMESHH
#define PZBUILDMULTIPHYSICSMESHH

#include <iostream>

#include "pzcompel.h"

/*
 *@brief This class has methods to build the mesh multiphysics
 */

class TPZBuildMultiphysicsMesh{
	
	
public:
	TPZBuildMultiphysicsMesh();
	
	~TPZBuildMultiphysicsMesh();
	
	/*
	 *@brief Creating multiphysic elements into mphysics computational mesh
	 Method to add elements in the mesh multiphysics
	 *@param cmeshVec [in]:  pointer to an vector of meshes
	 @param MFMesh [out]: my mesh multiphysics  
	 */
	void AddElements(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh);
	
	/*
	 *@brief Method to add connects in the mesh multiphisics
	 *@param cmeshVec [in]: pointer to an vector of meshes
	 @param MFMesh [out]: my mesh multiphysics  
	 */
	void AddConnects(TPZVec<TPZCompMesh *> cmeshVec, TPZCompMesh *MFMesh);
	
	/**
	 * @brief Transfer information from a specific set of meshes for the current mesh multiphysics
	 * @param meshVec[in] vector of meshes. Transfers the information
	 * @param MFMesh [out] mesh pointer who will receive the information
	 */	
	void TransferFromMeshes(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);
	
	/**
	 * @brief Transfer information from a specific mesh multiphysics for the current specific set of meshes 
	 * @param meshVec[out] vector of meshes that will receive the information.
	 * @param MFMesh [in] mesh pointer that Transfers the information 
	 */	
	void TransferFromMultiPhysics(TPZVec<TPZCompMesh *> &cmeshVec, TPZCompMesh *MFMesh);
	
	
};


#endif