cmake_minimum_required (VERSION 3.5)
project(pyClustering)

set(CMAKE_PREFIX_PATH "$ENV{VIRTUAL_ENV}")


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

include_directories(${CMAKE_PREFIX_PATH}/lib64/python3.12/site-packages/pybind11/include/pybind11)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(/usr/lib64/)

set(pybind11_DIR ${CMAKE_PREFIX_PATH}/lib64/python3.12/site-packages/pybind11/share/cmake/pybind11)
find_package(pybind11 REQUIRED)

include(cmake/CheckGit.cmake)
CheckGitSetup()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB MAC_FILES 
${PROJECT_SOURCE_DIR}/runSingle.py
${PROJECT_SOURCE_DIR}/sugarPosTest.csv
${PROJECT_SOURCE_DIR}/test.csv
${PROJECT_SOURCE_DIR}/testFile.root
${PROJECT_SOURCE_DIR}/__init__.py
)

foreach(_script ${MAC_FILES})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()

pybind11_add_module(clustering pyClustering.cc  ${sources} ${headers})
