cmake_minimum_required(VERSION 3.25)
set (CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-writable-strings")
project(QDD)


get_filename_component(USER $ENV{HOME} NAME)



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


add_subdirectory(src)
add_subdirectory(test)
