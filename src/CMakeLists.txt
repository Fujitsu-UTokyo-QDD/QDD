add_library(engine STATIC dd.cpp)
target_include_directories(engine PUBLIC ${PROJECT_SOURCE_DIR}/include ${Boost_INCLUDE_DIR})
target_include_directories(engine PUBLIC ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0 ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0/unsupported)
if(isMT)
  target_link_libraries(engine PUBLIC TBB::tbb Threads::Threads)
endif()
if(isMPI OR isMT)
  target_link_libraries(engine PUBLIC ${Boost_LIBRARIES})
endif()

add_library(alg STATIC algorithms/grover.cpp algorithms/shor.cpp)
target_include_directories(alg PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(alg PUBLIC engine task)


add_library(task STATIC task.cpp)
target_include_directories(task PUBLIC  ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(task PUBLIC  engine)
