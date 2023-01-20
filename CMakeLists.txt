cmake_minimum_required(VERSION 3.25)
set (CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-writable-strings")
project(QDD)

option(with_graphviz      "Use graphviz to draw task graph for debugging" OFF)


if(with_graphviz)
    find_package(PkgConfig REQUIRED)
    set(ENV{PKG_CONFIG_PATH} /Users/lsw/Code/graphviz/build/install/lib/pkgconfig)
    pkg_check_modules(CGRAPH REQUIRED IMPORTED_TARGET libcgraph)
    pkg_check_modules(GVC REQUIRED IMPORTED_TARGET libgvc)
endif()

if(CGRAPH_FOUND)
    message(STATUS ${CGRAPH_INCLUDEDIR})
    message(STATUS ${CGRAPH_LIBDIR})
    add_compile_definitions(GRAPHVIZ)
endif()


find_package(TBB REQUIRED)
add_library(engine STATIC dd.cpp  graph.cpp)
target_link_libraries(engine PUBLIC TBB::tbb )



add_executable(qdd main.cpp)

target_link_libraries(qdd PUBLIC engine)


if(with_graphviz)
    target_link_libraries(engine PUBLIC PkgConfig::CGRAPH PkgConfig::GVC)
endif()

