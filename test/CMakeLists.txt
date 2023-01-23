
add_executable(qdd_test test.cpp)

target_link_libraries(qdd_test PUBLIC engine PUBLIC GTest::gtest_main)

set(EIGEN_PATH ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0)

if(NOT EXISTS ${EIGEN_PATH})
    file(DOWNLOAD https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip "${EIGEN}.zip")
    file(ARCHIVE_EXTRACT INPUT "${EIGEN}.zip" DESTINATION ${PROJECT_SOURCE_DIR}/lib)
endif()



target_include_directories(qdd_test PRIVATE ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0 ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0/unsupported)


include(GoogleTest)
gtest_discover_tests(qdd_test)

add_executable(benchmark main.cpp)
target_link_libraries(benchmark PUBLIC engine alg)

