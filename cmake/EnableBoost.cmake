function(enable_boost target)
    find_package( Boost COMPONENTS graph unit_test_framework date_time  REQUIRED )
    target_link_libraries(${target} PRIVATE Boost::date_time Boost::graph)
    target_include_directories(${target} PRIVATE ${Boost_INCLUDE_DIRS})
    target_compile_definitions(${target} PRIVATE USING_BOOST)
    target_compile_definitions(${target} INTERFACE PZ_USING_BOOST)
    mark_as_advanced(Boost_GRAPH_LIBRARY_RELEASE  Boost_DATE_TIME_LIBRARY_RELEASE Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE Boost_INCLUDE_DIRS Boost_INCLUDE_DIR)
endfunction()
