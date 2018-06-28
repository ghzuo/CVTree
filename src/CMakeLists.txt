cmake_minimum_required (VERSION 2.6)
project (CVTree)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

ADD_DEFINITIONS("-O3 -std=c++11")

INCLUDE_DIRECTORIES("/usr/local/include")
LINK_DIRECTORIES(${LIBRARY_OUTPUT_PATH} "/usr/local/lib")
LINK_LIBRARIES(libcvtree.a z netcdf netcdf_c++)

SET(HEADFILES  global.h
  readgenome.h kstring.h
  stringOpt.h  memory.h
  neighborJoint.h distmatrix.h)

SET(LIBCVTREE_SRC ${HEADFILES}
  readgenome.cpp kstring.cpp
  stringOpt.cpp  memory.cpp
  neighborJoint.cpp distmatrix.cpp
  tree.cpp )
  
SET(CV_SRC ${HEADFILES}
  cvtree.h cvtree.cpp )

SET(CVDUMP_SRC ${HEADFILES}
    cvdump.cpp )

SET(TREE_SRC ${HEADFILES}
  buildTree.h buildTree.cpp)

ADD_LIBRARY(cvtree STATIC ${LIBCVTREE_SRC})

ADD_EXECUTABLE(cv ${CV_SRC})
ADD_DEPENDENCIES(cv cvtree)

ADD_EXECUTABLE(cvdump ${CVDUMP_SRC})
ADD_DEPENDENCIES(cvdump cvtree)

ADD_EXECUTABLE(tree ${TREE_SRC})
ADD_DEPENDENCIES(tree cvtree)

# INSTALL(TARGETS cvtree 
#   ARCHIVE DESTINATION lib)
# INSTALL(FILES ${HEADFILES}
#   DESTINATION include/cvtree)

INSTALL(TARGETS cv tree cvdump
  RUNTIME DESTINATION bin)
