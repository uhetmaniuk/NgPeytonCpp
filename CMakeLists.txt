PROJECT(CPPSELINV CXX)
cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

#--- Specify the include directories
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src"
                     "${PROJECT_SOURCE_DIR}/examples")

ADD_EXECUTABLE(test examples/clap2d5pt.cpp src/SymmetricSparse.h src/SymmetricSparse_impl.h
        examples/nd2d.h)

