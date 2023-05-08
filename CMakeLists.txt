cmake_minimum_required(VERSION 3.20)

project(ringsnark LANGUAGES CXX)

# Executables
add_executable(ringsnark-example)
target_sources(ringsnark-example PRIVATE ${CMAKE_CURRENT_LIST_DIR}/example.cpp)
target_link_libraries(ringsnark-example PRIVATE ringsnark)

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

find_package(Boost 1.80 COMPONENTS math)

## lib: polytools
target_link_libraries(ringsnark INTERFACE polytools)
target_include_directories(ringsnark INTERFACE ${Boost_INCLUDE_DIRS})
target_link_libraries(ringsnark INTERFACE ${Boost_LIBRARIES})
