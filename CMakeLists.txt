cmake_minimum_required(VERSION 3.6)
project(RELI)

set(CMAKE_CXX_STANDARD 11)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")

find_package(GSL REQUIRED)

include_directories (include)

set(SOURCE_FILES
        include/RELI_impl.h
        src/RELI_impl.cpp
        include/parser.h
        src/parser.cpp
        include/cxxopts.hpp
        include/utils.h
        src/RELI.cpp
        include/analysis.h)

add_executable(RELI ${SOURCE_FILES})
target_link_libraries(RELI GSL::gsl)
