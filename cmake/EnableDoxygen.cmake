function(enable_doxygen target)
    if (NOT "$ENV{DOXYGEN_ROOT}" STREQUAL "")
        message(STATUS "Looking for doxygen in $ENV{BOOST_ROOT}")
    endif()
    find_package(Doxygen
        REQUIRED dot
        REQUIRED)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs_doxygen/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.out)
    set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs_doxygen)
    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT})

    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
endfunction()
