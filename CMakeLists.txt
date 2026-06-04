cmake_minimum_required(VERSION 3.25)
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C CXX)


set(EIGEN_PATH ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0)

if(NOT EXISTS ${EIGEN_PATH})
    file(DOWNLOAD https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip "${EIGEN}.zip")
    file(ARCHIVE_EXTRACT INPUT "${EIGEN}.zip" DESTINATION ${PROJECT_SOURCE_DIR}/lib)
endif()
add_compile_definitions(EIGEN_DONT_VECTORIZE)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_compile_definitions(_GNU_SOURCE)
endif()
if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND NOT MSVC)
        add_compile_options(-fno-omit-frame-pointer)
    endif()
endif()

# clang-cl (LLVM on Windows with MSVC ABI)
if(MSVC AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_definitions(_USE_MATH_DEFINES)
    add_compile_options(/O2 /fp:fast /clang:-flto=thin)
    add_link_options(/clang:-flto=thin)
endif()

enable_testing()

include(FetchContent)
# On Windows with MSVC ABI (including clang-cl), GoogleTest must link against
# the shared CRT to match the runtime used by the rest of the build.
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.13.0
  FIND_PACKAGE_ARGS NAMES GTest
)
FetchContent_MakeAvailable(googletest)
# Mark GoogleTest headers as SYSTEM to suppress warnings from its headers
foreach(_gtest_target gtest gtest_main gmock gmock_main)
  if(TARGET ${_gtest_target})
    get_target_property(_inc ${_gtest_target} INTERFACE_INCLUDE_DIRECTORIES)
    if(_inc)
      set_target_properties(${_gtest_target} PROPERTIES
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_inc}")
    endif()
  endif()
endforeach()

option(isMPI "Configure for using MPI" OFF)
if (isMPI)
  add_compile_definitions(isMPI)
  find_package(MPI REQUIRED COMPONENTS CXX C)
  find_package(Boost REQUIRED COMPONENTS serialization mpi)
  message(STATUS "Boost include: ${Boost_INCLUDE_DIRS}")
endif()

add_subdirectory(src)
add_subdirectory(test)

option(BINDINGS "Configure for building Python bindings" ON)
if (BINDINGS)
  find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE nanobind_ROOT
  )
  find_package(nanobind CONFIG REQUIRED)
  add_subdirectory(qdd)
endif ()
