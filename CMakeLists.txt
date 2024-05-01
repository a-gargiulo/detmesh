cmake_minimum_required(VERSION 3.28)

project(detmesh VERSION 1.0.0)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED TRUE)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/build")

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Clang detected. Exporting compile commands...")
    set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
    add_compile_options(-Wall -Wextra -Weverything -pedantic)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/Wall)
endif()

file(GLOB SRCS "${PROJECT_SOURCE_DIR}/src/*.c")

add_executable(detmesh ${SRCS})

add_subdirectory("${PROJECT_SOURCE_DIR}/lib/fsolve")

target_include_directories(detmesh
    PUBLIC "${PROJECT_SOURCE_DIR}/lib/gmsh/include/"
    PUBLIC "${PROJECT_SOURCE_DIR}/lib/fsolve/include/"
    PUBLIC "${PROJECT_SOURCE_DIR}/include/"
)

target_link_directories(detmesh
    PUBLIC "${PROJECT_SOURCE_DIR}/lib/gmsh/"
)

target_link_libraries(detmesh
    PUBLIC gmsh
    PUBLIC fsolve
)

