cmake_minimum_required(VERSION 3.20)

project(ringsnark LANGUAGES CXX)


add_executable(bench)
target_sources(bench PRIVATE ${CMAKE_CURRENT_LIST_DIR}/microbench.cpp)
target_link_libraries(bench PRIVATE ringsnark)

# Includes
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/depends")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ringsnark")

target_include_directories(ringsnark INTERFACE ${CMAKE_CURRENT_LIST_DIR}/include/)

include_directories(.)

# Libraries
## lib: SEAL
# find_package(SEAL 4.0 REQUIRED)
if (TARGET SEAL::seal)
    target_link_libraries(ringsnark INTERFACE SEAL::seal)
elseif (TARGET SEAL::seal_shared)
    target_link_libraries(ringsnark INTERFACE SEAL::seal_shared)
else ()
    message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
endif ()

## lib: polytools
target_link_libraries(ringsnark INTERFACE polytools)

## lib: Boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.80 COMPONENTS math)

target_include_directories(ringsnark INTERFACE ${Boost_INCLUDE_DIRS})
target_link_libraries(ringsnark INTERFACE ${Boost_LIBRARIES})

## lib: Google benchmark
target_link_libraries(ringsnark INTERFACE benchmark::benchmark)

## lib: OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP enabled")
    target_link_libraries(ringsnark INTERFACE OpenMP::OpenMP_CXX)
endif()

# Features
include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if(result)
    message(STATUS "IPO / LTO enabled")
    set_target_properties(foo PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO enabled")
endif()

# Executables
## Examples
file( GLOB_RECURSE EXAMPLES_SOURCES RELATIVE "${CMAKE_SOURCE_DIR}/examples/" "examples/*.cpp" )
foreach( source_file ${EXAMPLES_SOURCES} )
    string( REPLACE ".cpp" "" source_name ${source_file} )
    add_executable( ${source_name} "examples/${source_file}" )
    target_link_libraries(${source_name} PRIVATE ringsnark)
endforeach( source_file ${EXAMPLES_SOURCES} )

## Benchmarks
file( GLOB_RECURSE EXAMPLES_SOURCES RELATIVE "${CMAKE_SOURCE_DIR}/benchmarks/" "benchmarks/*.cpp" )
foreach( source_file ${EXAMPLES_SOURCES} )
    string( REPLACE ".cpp" "" source_name ${source_file} )
    add_executable( ${source_name} "benchmarks/${source_file}" )
    target_link_libraries(${source_name} PRIVATE ringsnark)
endforeach( source_file ${EXAMPLES_SOURCES} )