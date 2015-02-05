SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

ADD_DEFINITIONS("-O3 -std=c++11")
LINK_LIBRARIES(z netcdf netcdf_c++)
  
SET(CV_SRC
  readgenome.h readgenome.cpp
  kstring.h kstring.cpp
  cvtree.h cvtree.cpp
  stringOpt.h stringOpt.cpp
  global.h
  )

SET(TREE_SRC
  buildTree.h buildTree.cpp
  stringOpt.h stringOpt.cpp
  memory.h memory.cpp
  distmatrix.h distmatrix.cpp
  kstring.h kstring.cpp
  tree.h tree.cpp
  neighborJoint.h neighborJoint.cpp
  global.h
  )


ADD_EXECUTABLE(cv ${CV_SRC})
ADD_EXECUTABLE(tree ${TREE_SRC})

INSTALL(TARGETS cv tree
  RUNTIME DESTINATION bin)
