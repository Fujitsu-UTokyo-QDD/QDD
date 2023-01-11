cmake_minimum_required(VERSION 3.25)
set (CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(QDD)


find_package(TBB REQUIRED)
add_library(engine STATIC dd.cpp engine.cpp complex.cpp epoch.cpp)
target_link_libraries(engine PUBLIC TBB::tbb )






add_executable(qdd main.cpp)

target_link_libraries(qdd PRIVATE engine)
