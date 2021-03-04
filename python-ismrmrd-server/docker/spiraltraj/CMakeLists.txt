cmake_minimum_required(VERSION 3.10)

# project name
PROJECT(spiraltraj)


FIND_PACKAGE(PythonInterp 3.6 REQUIRED)
FIND_PACKAGE(PythonLibs 3.6 REQUIRED)
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
include_directories(${PYTHON_INCLUDE_DIRS})

# This comes to hand if we also need to use the NumPy C API
exec_program(${PYTHON_EXECUTABLE}
             ARGS "-c \"import numpy; print(numpy.get_include())\""
             OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
             RETURN_VALUE NUMPY_NOT_FOUND
            )
if(NUMPY_NOT_FOUND)
    message(FATAL_ERROR "NumPy headers not found")
endif()

# compiler flags
#SET(CMAKE_CXX_FLAGS "-Wall -g -fPIC")
#MESSAGE("-- CXX Flags: " ${CMAKE_CXX_FLAGS})


set(SRCS "nonCartesianTraj.cpp" "vdspiral.cpp" "gen_spiral.cpp")
add_library(spiraltraj SHARED ${SRCS})


set_target_properties(spiraltraj
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME "spiraltraj"
        LINKER_LANGUAGE CXX
    )

#PYTHON_ADD_MODULE(gen_spiral gen_spiral.cpp)
#PYTHON_ADD_MODULE(gen_spiral ${SRCS})

#FILE(COPY hello.py DESTINATION .)
#ADD_TEST(NAME 01-HelloWorld COMMAND ${PYTHON_EXECUTABLE} hello.py)