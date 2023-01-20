
find_package(TBB REQUIRED)
add_library(engine STATIC dd.cpp  graph.cpp)
target_include_directories(engine PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(engine PUBLIC TBB::tbb)

if(with_graphviz)
    target_link_libraries(engine PUBLIC PkgConfig::CGRAPH PkgConfig::GVC)
endif()
