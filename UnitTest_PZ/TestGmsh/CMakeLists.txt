# @file neopz/UnitTest_PZ/TestGmsh/CMakeLists.txt  -- CMake file for unit test of reading gmsh files
configure_file("multiple_periodic_2d.msh" "multiple_periodic_2d.msh" COPYONLY)
configure_file("periodic_layer.msh" "periodic_layer.msh" COPYONLY)
add_unit_test(TestGmsh TestGmsh.cpp)
