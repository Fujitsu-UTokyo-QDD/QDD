add_executable(qdd_test test.cpp test_performance.cpp)
target_link_libraries(qdd_test PUBLIC engine PUBLIC GTest::gtest_main)
target_include_directories(qdd_test PRIVATE ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0 ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0/unsupported)
include(GoogleTest)
gtest_discover_tests(qdd_test)

if(isMPI)
  add_executable(mpi_test mpi_test.cpp)
  target_link_libraries(mpi_test PUBLIC engine PUBLIC GTest::gtest_main)
  gtest_discover_tests(mpi_test)

  add_executable(mpi_test_grover mpi_test_grover.cpp)
  target_link_libraries(mpi_test_grover PUBLIC engine)

  add_executable(serialization_test serialization_test.cpp)
  target_link_libraries(serialization_test PUBLIC engine)
  target_link_libraries(serialization_test PUBLIC engine PUBLIC GTest::gtest_main)
  gtest_discover_tests(serialization_test)
endif()
