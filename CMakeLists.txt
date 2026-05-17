cmake_minimum_required(VERSION 4.0)
set(CMAKE_INSTALL_PREFIX ../install/ CACHE PATH "" FORCE)
project(opti_tools LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -march=native")

add_library(optiToolsLib STATIC)
target_include_directories(optiToolsLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib ${CMAKE_CURRENT_SOURCE_DIR}/extern)
target_sources(optiToolsLib PRIVATE lib/solver/euler.cpp)

include(FetchContent)
FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG 12.1.0)
FetchContent_MakeAvailable(fmt)
target_link_libraries(optiToolsLib fmt::fmt)


if(${BUILD_FOR} MATCHES "PYTHON")
        find_package(Python 3.12 COMPONENTS Interpreter Development.Module REQUIRED)

        add_subdirectory(extern/nanobind)
        nanobind_add_module(opti_tools src/module/module.cpp)

        target_link_libraries(opti_tools PRIVATE optiToolsLib)

        set_target_properties(opti_tools PROPERTIES
                CXX_VISIBILITY_PRESET hidden
                VISIBILITY_INLINES_HIDDEN ON
        )

        nanobind_add_stub(
                opti_tools_stub
                MODULE opti_tools
                OUTPUT opti_tools.pyi
                PYTHON_PATH $<TARGET_FILE_DIR:opti_tools>
        )

        install(TARGETS opti_tools DESTINATION .)
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/opti_tools.pyi DESTINATION .)
else()
        add_executable(cpp_tests src/ascent_profile_drag.cpp)
        target_link_libraries(cpp_tests PUBLIC optiToolsLib)
        add_subdirectory(extern/matplotplusplus)
        target_link_libraries(cpp_tests PUBLIC matplot)
endif()
