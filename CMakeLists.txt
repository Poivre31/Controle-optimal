cmake_minimum_required(VERSION 4.0)
set(CMAKE_INSTALL_PREFIX ../install/ CACHE PATH "" FORCE)
project(controle_optimal LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -march=native")

add_library(math_lib STATIC)
target_include_directories(math_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib ${CMAKE_CURRENT_SOURCE_DIR}/extern)
target_sources(math_lib PRIVATE lib/lib.cpp)

include(FetchContent)
FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG 12.1.0)
FetchContent_MakeAvailable(fmt)
target_link_libraries(math_lib fmt::fmt)


add_executable(cpp_tests src/ascent_profile_drag.cpp)
target_link_libraries(cpp_tests PUBLIC math_lib)
add_subdirectory(extern/matplotplusplus)
target_link_libraries(cpp_tests PUBLIC matplot)
