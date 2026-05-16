cmake_minimum_required(VERSION 3.25)
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C CXX)


get_filename_component(USER $ENV{HOME} NAME)


set(EIGEN_PATH ${PROJECT_SOURCE_DIR}/lib/eigen-3.4.0)

if(NOT EXISTS ${EIGEN_PATH})
    file(DOWNLOAD https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip "${EIGEN}.zip")
    file(ARCHIVE_EXTRACT INPUT "${EIGEN}.zip" DESTINATION ${PROJECT_SOURCE_DIR}/lib)
endif()
add_compile_definitions(EIGEN_DONT_VECTORIZE)

add_compile_definitions(_GNU_SOURCE)
if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(-fno-omit-frame-pointer)
    endif()
endif()

enable_testing()

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.13.0
  FIND_PACKAGE_ARGS NAMES GTest
)
FetchContent_MakeAvailable(googletest)

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
  find_package(Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE nanobind_ROOT
  )
  find_package(nanobind CONFIG REQUIRED)
  add_subdirectory(qdd)
endif ()
