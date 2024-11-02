# @file neopz/Material/CMakeLists.txt  -- CMake file for the Material module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Material>
                           )

set(public_headers
    TPZMaterial.h
    TPZMaterialT.h
    TPZMatBase.h
    TPZBndCond.h
    TPZBndCondT.h
    TPZBndCondBase.h
    TPZMaterialData.h
    TPZMaterialDataT.h
    TPZMatTypes.h
    #generalinterfaces
    TPZMatError.h
    TPZMatWithMem.h
    TPZMatLoadCases.h
    TPZMatGeneralisedEigenVal.h
    TPZMatQuadraticEigenVal.h
    #single space materials
    TPZMatSingleSpace.h
    TPZMatErrorSingleSpace.h
    TPZMatInterfaceSingleSpace.h
    TPZMatTransientSingleSpace.h
    TPZMatSingleSpaceBCSol.h
    #combined spaces materials
    TPZMatCombinedSpaces.h
    TPZMatErrorCombinedSpaces.h
    TPZMatInterfaceCombinedSpaces.h
    #other
    TPZNullMaterial.h
    TPZNullMaterialCS.h
    TPZLagrangeMultiplier.h
    TPZLagrangeMultiplierCS.h
    TPZMatTypes.h
    TPZMatInterfaceHybridElasticityStokes.h
    )

set(headers
    TPZMaterial.h
    TPZMaterialT.h
    TPZMatBase.h
    TPZBndCond.h
    TPZBndCondT.h
    TPZBndCondBase.h
    TPZMaterialData.h
    TPZMaterialDataT.h
    TPZMatTypes.h
    #generalinterfaces
    TPZMatError.h
    TPZMatWithMem.h
    TPZMatLoadCases.h
    TPZMatGeneralisedEigenVal.h
    TPZMatQuadraticEigenVal.h
    #single space materials
    TPZMatSingleSpace.h
    TPZMatErrorSingleSpace.h
    TPZMatInterfaceSingleSpace.h
    TPZMatTransientSingleSpace.h
    TPZMatSingleSpaceBCSol.h
    #combined spaces materials
    TPZMatCombinedSpaces.h
    TPZMatErrorCombinedSpaces.h
    TPZMatInterfaceCombinedSpaces.h
    #other
    TPZNullMaterial.h
    TPZNullMaterialCS.h
    TPZLagrangeMultiplier.h
    TPZLagrangeMultiplierCS.h
    TPZMatInterfaceHybridElasticityStokes.h
   )
set(sources
    TPZMaterial.cpp
    TPZMaterialT.cpp
    TPZBndCond.cpp
    TPZBndCondT.cpp
    TPZMaterialData.cpp
    TPZMaterialDataT.cpp
    #general interfaces
    TPZMatError.cpp
    TPZMatLoadCases.cpp
    TPZMatGeneralisedEigenVal.cpp
    TPZMatQuadraticEigenVal.cpp
    #single space materials
    TPZMatSingleSpace.cpp
    TPZMatErrorSingleSpace.cpp
    TPZMatInterfaceSingleSpace.cpp
    TPZMatTransientSingleSpace.cpp
    TPZMatSingleSpaceBCSol.cpp
    #combined spaces materials
    TPZMatCombinedSpaces.cpp
    TPZMatErrorCombinedSpaces.cpp
    TPZMatInterfaceCombinedSpaces.cpp
    #other
    TPZNullMaterial.cpp
    TPZNullMaterialCS.cpp
    TPZLagrangeMultiplier.cpp
    TPZLagrangeMultiplierCS.cpp
    TPZMatInterfaceHybridElasticityStokes.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

if(BUILD_PLASTICITY_MATERIALS)
    if ((NOT STATE_TYPE STREQUAL "double") OR (NOT  REAL_TYPE STREQUAL "double"))
		message (FATAL_ERROR
            "\nConflicting options for BUILD_PLASTICITY_MATERIALS\n"
            "Set REAL_TYPE=double and STATE_TYPE=double to build them\n"
            "Or set BUILD_PLASTICITY_MATERIALS=OFF")
    else()
        add_subdirectory(Plasticity)
    endif()
endif()
#elasticity materials
add_subdirectory(Elasticity)
#materials for conservation laws
add_subdirectory(ConsLaw)
#materials for black oil
add_subdirectory(BlackOil)
#projection materials
add_subdirectory(Projection)
#poisson materials
add_subdirectory(Poisson)
#electromagnetic materials
add_subdirectory(Electromagnetics)
#Darcy flow materials
add_subdirectory(DarcyFlow)

set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)