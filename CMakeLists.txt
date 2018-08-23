###
# Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai, China.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2017-09-01 16:22:53
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2018-04-26 19:50:44
###

cmake_minimum_required(VERSION 2.6)
project (CVTree)

## For Compile
##########################################################
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(CMAKE_VERBOSE_MAKEFILE on)

## use C++11
ADD_DEFINITIONS("-O3 -std=c++11")

### Check the OPENMP
find_package(OpenMP)
if (OPENMP_FOUND)
   set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

### Check require LIABRARIES
INCLUDE_DIRECTORIES("/usr/local/include")
LINK_DIRECTORIES(${LIBRARY_OUTPUT_PATH} "/usr/local/lib")
set(REQUIRED_LIBRARIES ${REQUIRED_LIBRARIES} z)
set(REQUIRED_LIBRARIES ${REQUIRED_LIBRARIES} netcdf_c++)
foreach(LIB ${REQUIRED_LIBRARIES})
  find_library(THE${LIB} ${LIB})
  if(THE${LIB})
    message("-- Found library " ${LIB} ": " ${THE${LIB}})
  else()
    message(FATAL_ERROR ${LIB} " library not found!")
  endif()
endforeach()
LINK_LIBRARIES(${REQUIRED_LIBRARIES})

## include subdirectories
ADD_SUBDIRECTORY(src)

## for the latex manual
find_package(LATEX)
if(LATEX_FOUND)
  ADD_SUBDIRECTORY(manual)
endif()



