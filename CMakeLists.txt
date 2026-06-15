cmake_minimum_required(VERSION 4.0)
set(CMAKE_INSTALL_PREFIX ../install/ CACHE PATH "" FORCE)
project(controle_optimal LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(math_lib STATIC)
target_include_directories(math_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib ${CMAKE_CURRENT_SOURCE_DIR}/extern)
target_sources(math_lib PRIVATE lib/lib.cpp)
target_compile_options(math_lib PRIVATE -O3 -march=native)

include(FetchContent)
fetchcontent_declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG 12.1.0)
fetchcontent_makeavailable(fmt)
target_link_libraries(math_lib fmt::fmt)


add_executable(cpp_tests src/ascent_profile_no_drag.cpp)
target_link_libraries(cpp_tests PUBLIC math_lib)
add_subdirectory(extern/matplotplusplus)
target_link_libraries(cpp_tests PUBLIC matplot)

target_compile_options(cpp_tests PRIVATE -O3 -march=native)
