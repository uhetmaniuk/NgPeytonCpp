cmake_minimum_required(VERSION 3.12)
project(NgPeytonCpp CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_options(-Wall -Wextra)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

include_directories("${PROJECT_SOURCE_DIR}/src"
                     "${PROJECT_SOURCE_DIR}/examples")

add_executable(clap2d5pt examples/clap2d5pt.cpp)
add_executable(slap2d9pt examples/slap2d9pt.cpp)
add_executable(dlap2d9pt examples/dlap2d9pt.cpp)

# --- Unit tests (ctest) ---
enable_testing()

add_executable(test_isFullRepresentation tests/test_isFullRepresentation.cpp)
target_include_directories(test_isFullRepresentation PRIVATE "${PROJECT_SOURCE_DIR}/src")
add_test(NAME isFullRepresentation COMMAND test_isFullRepresentation)

add_executable(test_solve tests/test_solve.cpp)
target_include_directories(test_solve PRIVATE "${PROJECT_SOURCE_DIR}/src")
add_test(NAME solve COMMAND test_solve)

add_test(NAME clap2d5pt COMMAND clap2d5pt)
