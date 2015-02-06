SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

ADD_DEFINITIONS("-O3 -std=c++11")

LINK_DIRECTORIES(${LIBRARY_OUTPUT_PATH})
LINK_LIBRARIES(z netcdf netcdf_c++ libcvtree.a)

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

SET(TREE_SRC ${HEADFILES}
  buildTree.h buildTree.cpp)

ADD_LIBRARY(cvtree STATIC ${LIBCVTREE_SRC})

ADD_EXECUTABLE(cv ${CV_SRC})
ADD_DEPENDENCIES(cv cvtree)

ADD_EXECUTABLE(tree ${TREE_SRC})
ADD_DEPENDENCIES(tree cvtree)

# INSTALL(TARGETS cvtree 
#   ARCHIVE DESTINATION lib)
# INSTALL(FILES ${HEADFILES}
#   DESTINATION include/cvtree)

INSTALL(TARGETS cv tree
  RUNTIME DESTINATION bin)
