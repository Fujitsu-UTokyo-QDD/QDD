cmake_minimum_required(VERSION 3.25)
set (CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(QDD)

add_library(CLQueue INTERFACE "chase_lev_queue.hpp")

target_compile_features(CLQueue INTERFACE cxx_std_20)

find_package(TBB REQUIRED)
add_library(engine STATIC dd.cpp engine.cpp  epoch.cpp)
target_link_libraries(engine PUBLIC TBB::tbb  CLQueue)






add_executable(qdd main.cpp)

target_link_libraries(qdd PRIVATE engine)
