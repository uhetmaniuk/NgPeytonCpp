PROJECT(CPPSELINV CXX)
cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

#--- Specify the include directories
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src"
                     "${PROJECT_SOURCE_DIR}/examples")

ADD_EXECUTABLE(clap2d5pt examples/clap2d5pt.cpp src/SymmetricSparse.h src/SymmetricSparse_impl.h
        examples/nd2d.h)
ADD_EXECUTABLE(slap2d9pt examples/slap2d9pt.cpp)
ADD_EXECUTABLE(dlap2d9pt examples/dlap2d9pt.cpp)

# --- Unit tests (ctest) ---
enable_testing()

ADD_EXECUTABLE(test_isFullRepresentation tests/test_isFullRepresentation.cpp)
target_include_directories(test_isFullRepresentation PRIVATE "${PROJECT_SOURCE_DIR}/src")
ADD_TEST(NAME isFullRepresentation COMMAND test_isFullRepresentation)

ADD_EXECUTABLE(test_solve tests/test_solve.cpp)
target_include_directories(test_solve PRIVATE "${PROJECT_SOURCE_DIR}/src")
ADD_TEST(NAME solve COMMAND test_solve)

