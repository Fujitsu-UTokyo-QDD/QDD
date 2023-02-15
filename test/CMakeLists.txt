
add_executable(qdd_test test.cpp)

target_link_libraries(qdd_test PUBLIC engine PUBLIC GTest::gtest_main)




target_include_directories(qdd_test PRIVATE ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0 ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0/unsupported)


include(GoogleTest)
gtest_discover_tests(qdd_test)

add_executable(benchmark main.cpp)
target_link_libraries(benchmark PUBLIC engine alg)

add_executable(gc_test  qv15.cpp)
target_link_libraries(gc_test PUBLIC engine alg)

add_executable(fiber_test  fiber.cpp)
target_link_libraries(fiber_test PRIVATE task)

