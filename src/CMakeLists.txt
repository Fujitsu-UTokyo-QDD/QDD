
find_package(TBB REQUIRED)
add_library(engine STATIC dd.cpp  graph.cpp epoch.cpp)
target_include_directories(engine PUBLIC ${PROJECT_SOURCE_DIR}/include ${Boost_INCLUDE_DIR})
target_include_directories(engine PUBLIC ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0 ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0/unsupported)
target_link_libraries(engine PUBLIC TBB::tbb Threads::Threads ${Boost_LIBRARIES})

if(with_graphviz)
    target_link_libraries(engine PUBLIC PkgConfig::CGRAPH PkgConfig::GVC)
endif()




add_library(task STATIC task.cpp)
target_include_directories(task PUBLIC  ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(task PUBLIC  engine)

add_library(alg STATIC algorithms/grover.cpp algorithms/shor.cpp)
target_include_directories(alg PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(alg PUBLIC engine task)
