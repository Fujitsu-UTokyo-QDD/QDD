add_library(engine STATIC dd.cpp)
target_include_directories(engine PUBLIC ${PROJECT_SOURCE_DIR}/include ${Boost_INCLUDE_DIR})
target_include_directories(engine PUBLIC ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0 ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0/unsupported)
target_link_libraries(engine PUBLIC numa)
if(isMT)
  target_link_libraries(engine PUBLIC TBB::tbb Threads::Threads)
endif()
if(isMPI OR isMT)
  target_link_libraries(engine PUBLIC ${Boost_LIBRARIES})
endif()

if(isMT)
add_library(alg STATIC algorithms/grover.cpp algorithms/shor.cpp)
target_include_directories(alg PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(alg PUBLIC engine task)


add_library(task STATIC task.cpp)
target_include_directories(task PUBLIC  ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(task PUBLIC  engine)
endif()
#set(CMAKE_MODULE_PATH
#    "${PROJECT_SOURCE_DIR}/CacheLib/opt/cachelib/lib/cmake"
#  ${CMAKE_MODULE_PATH})
set(wangle_DIR "${PROJECT_SOURCE_DIR}/CacheLib/opt/cachelib/lib/cmake/wangle/")
set(mvfst_DIR "${PROJECT_SOURCE_DIR}/CacheLib/opt/cachelib/lib/cmake/mvfst/")
set(Fizz_DIR "${PROJECT_SOURCE_DIR}/CacheLib/opt/cachelib/lib/cmake/fizz/")
set(fizz_DIR "${PROJECT_SOURCE_DIR}/CacheLib/opt/cachelib/lib/cmake/fizz/")
set(FBThrift_DIR "${PROJECT_SOURCE_DIR}/CacheLib/opt/cachelib/lib/cmake/fbthrift/")
set(folly_DIR "${PROJECT_SOURCE_DIR}/CacheLib/opt/cachelib/lib/cmake/folly/")
set(fmt_DIR "${PROJECT_SOURCE_DIR}/CacheLib/opt/cachelib/lib/cmake/fmt/")
find_package(cachelib PATHS "${PROJECT_SOURCE_DIR}/CacheLib/opt/cachelib/lib/cmake/cachelib")
target_link_libraries(engine PUBLIC cachelib)
message(STATUS ${CACHELIB_INCLUDE_DIR})
target_include_directories(engine PUBLIC ${CACHELIB_INCLUDE_DIR})

if(CACHELIB)
    message(STATUS "Use cachelib opcache")
    target_compile_definitions(engine PRIVATE CACHELIB_OPCACHE)
endif()



