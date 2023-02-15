cmake_minimum_required(VERSION 3.25)
set (CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(QDD)


get_filename_component(USER $ENV{HOME} NAME)


set(EIGEN_PATH ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0)

if(NOT EXISTS ${EIGEN_PATH})
    file(DOWNLOAD https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip "${EIGEN}.zip")
    file(ARCHIVE_EXTRACT INPUT "${EIGEN}.zip" DESTINATION ${PROJECT_SOURCE_DIR}/lib)
endif()

option(with_mt         "Run with multiple threads" OFF)
if(with_mt)
    message(STATUS "Run with multiple threads")
    add_compile_definitions(MT)
endif()


option(with_graphviz      "Use graphviz to draw task graph for debugging" OFF)
if(with_graphviz)
    find_package(PkgConfig REQUIRED)
    if(${USER} STREQUAL "lsw")
        set(ENV{PKG_CONFIG_PATH} /Users/lsw/Code/graphviz/build/install/lib/pkgconfig)
    endif()
    pkg_check_modules(CGRAPH REQUIRED IMPORTED_TARGET libcgraph)
    pkg_check_modules(GVC REQUIRED IMPORTED_TARGET libgvc)
endif()

if(CGRAPH_FOUND)
    message(STATUS ${CGRAPH_INCLUDEDIR})
    message(STATUS ${CGRAPH_LIBDIR})
    add_compile_definitions(GRAPHVIZ)
endif()



set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
add_compile_definitions(_GUN_SOURCE)

enable_testing()

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.13.0
  FIND_PACKAGE_ARGS NAMES GTest
)
FetchContent_MakeAvailable(googletest)

find_package(Boost REQUIRED context)
message(STATUS "Boost include: ${Boost_INCLUDE_DIRS}")


add_subdirectory(src)
add_subdirectory(test)
