# @file neopz/Mesh/CMakeLists.txt  -- CMake file for the Mesh module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Mesh>
                           )


set(headers
  doxmesh.h
  pzconnect.h
  pzgeoelbc.h  
  TPZCompMeshTools.h
  TPZGeoMeshTools.h
  pzcheckgeom.h
  pzgeoel.h
  pzintel.h  
  TPZGeoElement.h
  TPZInterfaceEl.h
  pzcheckmesh.h
  pzelchdivbound2.h
  pzgeoelrefless.h
  pzinterpolationspace.h
  TPZCompElDisc.h
  TPZGeoElement.h.h
  TPZMultiphysicsCompMesh.h
  pzcheckrestraint.h
  pzelchdiv.h
  pzgeoelrefless.h.h
  pzmultiphysicscompel.h
  tpzgeoelmapped.h
  TPZMultiphysicsInterfaceEl.h
  pzcmesh.h
  pzelctemp.h
  pzgeoelside.h
  pzmultiphysicselement.h
  TPZCompElHCurl.h
  tpzgeoelrefpattern.h
  TPZOneShapeRestraint.h
  pzcompel.h
  pzelementgroup.h
  pzgmesh.h
  pzreducedspace.h
  TPZCompElHDivCollapsed.h
  tpzgeoelrefpattern.h.h
  pzcompelwithmem.h
  pzelmat.h
  TPZElementMatrixT.h
  pzgnode.h
  TPZCompElLagrange.h
  TPZGeoElSideAncestors.h
  pzcondensedcompel.h
  pzsubcmesh.h
  TPZGeoElSidePartition.h
  pzflowcmesh.h
  tpzagglomeratemesh.h
  TPZAgglomerateEl.h
  TPZCompElH1.h
  TPZCompElKernelHDiv.h
  TPZCompElKernelHDiv3D.h
  TPZHCurlEquationFilter.h
  TPZCompElDiscScaled.h
  TPZCompElHDivDuplConnects.h
  TPZCompElHDivDuplConnectsBound.h
  TPZCompElUnitaryLagrange.h
  SPZPeriodicData.h
)

set(sources
  pzcheckgeom.cpp
  pzelchdivbound2.cpp
  pzgeoelside.cpp
  pzreducedspace.cpp
  TPZGeoElSidePartition.cpp
  pzcheckmesh.cpp
  pzelchdiv.cpp
  pzgmesh.cpp
  TPZCompElHDivCollapsed.cpp
  TPZGeoMeshTools.cpp
  pzcheckrestraint.cpp
  pzelctemp.cpp
  pzgnode.cpp
  TPZCompElLagrange.cpp
  TPZInterfaceEl.cpp
  pzcmesh.cpp
  pzelementgroup.cpp  
  pzsubcmesh.cpp
  TPZMultiphysicsCompMesh.cpp
  pzcompel.cpp
  pzelmat.cpp
  TPZElementMatrixT.cpp
  TPZCompMeshTools.cpp
  TPZMultiphysicsInterfaceEl.cpp
  pzcompelwithmem.cpp  
  pzintel.cpp  
  TPZGeoElement.cpp
  TPZOneShapeRestraint.cpp
  pzcondensedcompel.cpp
  pzgeoelbc.cpp
  pzinterpolationspace.cpp
  TPZAgglomerateMesh.cpp
  tpzgeoelmapped.cpp
  pzconnect.cpp
  pzgeoel.cpp
  pzmultiphysicscompel.cpp
  TPZCompElDisc.cpp
  tpzgeoelrefpattern.cpp
  pzgeoelrefless.cpp
  pzmultiphysicselement.cpp
  TPZCompElHCurl.cpp
  TPZGeoElSideAncestors.cpp
  pzflowcmesh.cpp
  TPZAgglomerateEl.cpp
  TPZCompElH1.cpp
  TPZCompElKernelHDiv.cpp
  TPZCompElKernelHDiv3D.cpp
  TPZHCurlEquationFilter.cpp
  TPZCompElDiscScaled.cpp
  TPZCompElHDivDuplConnects.cpp
  TPZCompElHDivDuplConnectsBound.cpp
  TPZCompElUnitaryLagrange.cpp
)

if(USING_LAPACK)
  list(APPEND headers
    TPZSBFemElementGroup.h
    TPZSBFemMultiphysicsElGroup.h
    TPZSBFemVolume.h
	  TPZSBFemVolumeHdiv.h
	  TPZSBFemVolumeL2.h
	  TPZSBFemVolumeMultiphysics.h
    TPZCompElHDivSBFem.h)
  list(APPEND sources
    TPZSBFemElementGroup.cpp
    TPZSBFemMultiphysicsElGroup.cpp
    TPZSBFemVolume.cpp
	  TPZSBFemVolumeHdiv.cpp
	  TPZSBFemVolumeL2.cpp
	  TPZSBFemVolumeMultiphysics.cpp
    TPZCompElHDivSBFem.cpp)
endif()

  
install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Mesh)

target_sources(pz PRIVATE ${headers} ${sources})


#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)
