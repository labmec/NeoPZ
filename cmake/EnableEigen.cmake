# Here we will search for eigen. If eigen is not found, we will download it.
function(enable_eigen target)
    #perhaps eigen was already downloaded when running cmake
    if(NOT eigen_POPULATED)
        #find_package(Eigen3 3.4 QUIET CONFIG)

        if(NOT EIGEN3_FOUND)
            # Couldn't load via target, so fall back to allowing module mode finding, which will pick up
            # tools/FindEigen3.cmake
            # find_package(Eigen3 3.4 QUIET)
        endif()
        if(NOT EIGEN3_FOUND)
            set(EIGEN3_VERSION_STRING "7dd3dda3daa218147557b33f8d05b3b023f05f7d")
            include(FetchContent)
            FetchContent_Declare(
                eigen
                GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
                GIT_TAG ${EIGEN3_VERSION_STRING})

            FetchContent_GetProperties(eigen)
            if(NOT eigen_POPULATED)
                message(STATUS "Downloading Eigen")
                FetchContent_Populate(eigen)
                set(DOWNLOADED_EIGEN TRUE)
            endif()

            set(EIGEN3_INCLUDE_DIR ${eigen_SOURCE_DIR})
            message("EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
            set(EIGEN3_FOUND TRUE)
        endif()
    endif()
    if(EIGEN3_FOUND)
        # if eigen was downloaded, the target was not created
        if(NOT TARGET Eigen3::Eigen)
            add_library(Eigen3::Eigen IMPORTED INTERFACE)
            target_link_libraries(pz PUBLIC ${Eigen3_LIBRARIES})
            include_directories(${EIGEN3_INCLUDE_DIR})
			# target_link_libraries(${target} PUBLIC Eigen3::Eigen)
			target_compile_definitions(${target} PRIVATE USING_EIGEN)
            target_compile_definitions(${target} INTERFACE PZ_USING_EIGEN)
            # target_include_directories(${target} PUBLIC ${EIGEN3_INCLUDE_DIR})
            # set_property(TARGET Eigen3::Eigen PROPERTY INTERFACE_INCLUDE_DIRECTORIES
            #    "${EIGEN3_INCLUDE_DIR}")
        endif()

        # Eigen 3.3.1+ cmake sets EIGEN3_VERSION_STRING (and hard codes the version when installed
        # rather than looking it up in the cmake script); older versions, and the
        # tools/FindEigen3.cmake, set EIGEN3_VERSION instead.
        if(NOT EIGEN3_VERSION AND EIGEN3_VERSION_STRING)
            set(EIGEN3_VERSION ${EIGEN3_VERSION_STRING})
        endif()
        message(STATUS "Eigen found at ${EIGEN3_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR "Could not satisfy dependency: Eigen")
    endif()
endfunction()