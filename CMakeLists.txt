cmake_minimum_required(VERSION 3.12)
project(NgPeytonCpp CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_options(-Wall -Wextra)

# --- Library ---
add_library(NgPeytonCpp src/LDLtSolver.cpp)
target_include_directories(NgPeytonCpp
  PUBLIC "${PROJECT_SOURCE_DIR}/include"
  PRIVATE "${PROJECT_SOURCE_DIR}/src")

# --- Examples ---
add_executable(clap2d5pt examples/clap2d5pt.cpp)
target_include_directories(clap2d5pt PRIVATE "${PROJECT_SOURCE_DIR}/examples")
target_link_libraries(clap2d5pt NgPeytonCpp)

add_executable(slap2d9pt examples/slap2d9pt.cpp)
target_include_directories(slap2d9pt PRIVATE "${PROJECT_SOURCE_DIR}/examples")
target_link_libraries(slap2d9pt NgPeytonCpp)

add_executable(dlap2d9pt examples/dlap2d9pt.cpp)
target_include_directories(dlap2d9pt PRIVATE "${PROJECT_SOURCE_DIR}/examples")
target_link_libraries(dlap2d9pt NgPeytonCpp)

# --- Unit tests (ctest) ---
enable_testing()

add_executable(test_isFullRepresentation tests/test_isFullRepresentation.cpp)
target_include_directories(test_isFullRepresentation
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
add_test(NAME isFullRepresentation COMMAND test_isFullRepresentation)

add_executable(test_solve tests/test_solve.cpp)
target_link_libraries(test_solve NgPeytonCpp)
add_test(NAME solve COMMAND test_solve)

add_test(NAME clap2d5pt COMMAND clap2d5pt)
add_test(NAME dlap2d9pt COMMAND dlap2d9pt)
add_test(NAME slap2d9pt COMMAND slap2d9pt)
