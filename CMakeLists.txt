cmake_minimum_required(VERSION 3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT})
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(ZLIB REQUIRED)

# upstream flags minus -march=native and -flto so wheels stay portable
set(SANITY_FLAGS -O3 -fno-math-errno)

# upstream sources except the CLI translation unit, which _sanity.cpp includes
add_library(sanity_upstream STATIC
  src/FitFrac.cpp
  src/helpers/io/ReadInputFiles.cpp
  src/helpers/io/DenseRowReader.cpp
  src/helpers/io/npy_writer.cpp
  src/helpers/math/FukushimaLambertW.cc
  src/helpers/math/Digamma_Trigamma.cpp
)
target_include_directories(sanity_upstream PUBLIC src)
target_compile_options(sanity_upstream PRIVATE ${SANITY_FLAGS})
target_link_libraries(sanity_upstream PUBLIC OpenMP::OpenMP_CXX ZLIB::ZLIB)

set(PYX ${CMAKE_CURRENT_SOURCE_DIR}/src/pysanity/_core.pyx)
set(PYX_CXX ${CMAKE_CURRENT_BINARY_DIR}/_core.cxx)
add_custom_command(
  OUTPUT ${PYX_CXX}
  COMMAND ${Python_EXECUTABLE} -m cython --cplus -3 ${PYX} -o ${PYX_CXX}
  DEPENDS ${PYX} src/pysanity/_sanity.h
  VERBATIM
)

if(SKBUILD_SABI_VERSION)
  python_add_library(_core MODULE ${PYX_CXX} src/pysanity/_sanity.cpp WITH_SOABI USE_SABI ${SKBUILD_SABI_VERSION})
else()
  python_add_library(_core MODULE ${PYX_CXX} src/pysanity/_sanity.cpp WITH_SOABI)
endif()
target_include_directories(_core PRIVATE src/pysanity)
target_compile_options(_core PRIVATE ${SANITY_FLAGS})
target_link_libraries(_core PRIVATE sanity_upstream)

install(TARGETS _core DESTINATION pysanity)
