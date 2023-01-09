cmake_minimum_required(VERSION 3.25)
set (CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(QDD)


add_library(engine STATIC dd.cpp engine.cpp complex.cpp)
target_include_directories(engine PUBLIC /Users/lsw/Code/oneTBB/build/install/include)






add_executable(qdd main.cpp)

target_link_libraries(qdd PRIVATE engine)
