cmake_minimum_required(VERSION 3.24.2)
project("MeshSubdivider")
set (CMAKE_CXX_STANDARD 17)

include_directories("${CMAKE_CURRENT_LIST_DIR}/Header")
file(GLOB_RECURSE CPPFiles CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/Source/*.cpp")
file(GLOB_RECURSE HFiles CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/Header/*.h")
add_executable(MeshSubdivider ${CPPFiles} ${HFiles})
