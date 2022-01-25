function(enable_mkl target)
  find_package(MKL REQUIRED)

  set(MKL_THREAD_MODEL_TYPES "Options are: seq, tbb or omp(default)"
    "\nseq: Sequential\ntbb:TBB\nOMP:OpenMP")
  if(NOT MKL_THREAD_MODEL AND TARGET mkl::mkl_intel_32bit_tbb_dyn)
    set(MKL_THREAD_MODEL "tbb" CACHE STRING "${MKL_THREAD_MODEL_TYPES}")
  else()
    set(MKL_THREAD_MODEL "omp" CACHE STRING "${MKL_THREAD_MODEL_TYPES}")
    endif()
  
  if(NOT MKL_THREAD_MODEL STREQUAL "tbb" AND
      NOT MKL_THREAD_MODEL STREQUAL "omp" AND
      NOT MKL_THREAD_MODEL STREQUAL "seq")
    message(FATAL_ERROR "MKL THREADING MODEL ERROR"
      "\n${MKL_THREAD_MODEL_TYPES}")
  endif()
  message(STATUS "Setting MKL threading model: ${MKL_THREAD_MODEL}")
  
  #do NOT change this lib unless you know what you are doing
  target_link_libraries(${target} PRIVATE mkl::mkl_intel_32bit_${MKL_THREAD_MODEL}_dyn)

  if(MKL_THREAD_MODEL STREQUAL "tbb")
    include(cmake/EnableTBB.cmake)
	  enable_tbb(pz)
  endif()
  #on our test machine it was needed to link directly with mkl core
  #perhaps on newer mkl installs this is not needed anymore?
  #if(APPLE)
 #     target_link_libraries(${target} PRIVATE ${_mkl_core_lib})
  #endif()
  target_include_directories(${target} PRIVATE ${MKL_INCLUDE_DIR})
  target_compile_definitions(${target} PRIVATE USING_MKL)
  target_compile_definitions(${target} INTERFACE PZ_USING_MKL)

  set(USING_LAPACK ON CACHE PATH "Whether the LAPACK library will be linked in" FORCE)
  #TODOWIN32: should we do something with mkl_rt on windows?
  mark_as_advanced(MKL_THREAD_MODEL)
endfunction()
