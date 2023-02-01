
find_package(TBB REQUIRED)
add_library(engine STATIC dd.cpp  graph.cpp)
target_include_directories(engine PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(engine PUBLIC TBB::tbb Threads::Threads)
add_library(engineSHARED SHARED dd.cpp  graph.cpp)
target_include_directories(engineSHARED PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(engineSHARED PUBLIC TBB::tbb Threads::Threads)

if(with_graphviz)
    target_link_libraries(engine PUBLIC PkgConfig::CGRAPH PkgConfig::GVC)
    target_link_libraries(engineSHARED PUBLIC PkgConfig::CGRAPH PkgConfig::GVC)
endif()


add_library(alg STATIC algorithms/grover.cpp algorithms/shor.cpp)
target_include_directories(alg PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(alg PUBLIC engine)

add_library(algSHARED SHARED algorithms/grover.cpp algorithms/shor.cpp)
target_include_directories(algSHARED PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(algSHARED PUBLIC engineSHARED)
