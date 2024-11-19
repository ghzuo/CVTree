###
# Copyright (c) 2018  Wenzhou Institute, University of Chinese Academy of Sciences.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2017-09-01 15:55:03
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2024-11-19 17:12:22
###

INCLUDE_DIRECTORIES("../kit")

SET(CVTREEHEADS 
  genome.h 
  kstring.h cvmeth.h
  distmatrix.h distmeth.h 
  tree.h treemeth.h 
  sampling.h sampleMeth.h
  )

SET(LIBCVTREE_SRC ${CVTREEHEADS}
  genome.cpp
  kstring.cpp cvmeth.cpp
  distmeth.cpp distmatrix.cpp 
  tree.cpp treemeth.cpp 
  sampling.cpp sampleMeth.cpp
  )

SET(CV_SRC ${CVTREEHEADS}
  g2cv.h g2cv.cpp)

SET(CVDUMP_SRC ${CVTREEHEADS}
  cvdump.cpp)

SET(DIST_SRC ${CVTREEHEADS}
  cv2dm.h cv2dm.cpp)

SET(GETDIST_SRC ${CVTREEHEADS}
  getdist.cpp)

SET(MCONV_SRC ${CVTREEHEADS}
  mconv.cpp)

SET(DIFFMTX_SRC ${CVTREEHEADS} diffMatrix.cpp)

SET(TREEDIST_SRC ${CVTREEHEADS}
    treedist.h treedist.cpp
    setsym.h
    marktree.cpp marktree.h
  )

SET(TREEBOOT_SRC ${CVTREEHEADS}
  treeboot.h treeboot.cpp
  setsym.h
  marktree.cpp marktree.h
)

SET(NJ_SRC ${CVTREEHEADS}
  dm2tree.h dm2tree.cpp)

SET(CVTREE_SRC ${CVTREEHEADS}
  marktree.h marktree.cpp
  setsym.h
  cvtree.h cvtree.cpp)

ADD_LIBRARY(cvtr STATIC ${LIBCVTREE_SRC})

ADD_EXECUTABLE(g2cv ${CV_SRC})
TARGET_LINK_LIBRARIES(g2cv cvtr kit)

ADD_EXECUTABLE(cvdump ${CVDUMP_SRC})
TARGET_LINK_LIBRARIES(cvdump cvtr kit)

ADD_EXECUTABLE(cv2dm ${DIST_SRC})
TARGET_LINK_LIBRARIES(cv2dm cvtr kit)

ADD_EXECUTABLE(getdist ${GETDIST_SRC})
TARGET_LINK_LIBRARIES(getdist cvtr kit)

ADD_EXECUTABLE(mconv ${MCONV_SRC})
TARGET_LINK_LIBRARIES(mconv cvtr kit)

ADD_EXECUTABLE(diffmtx ${DIFFMTX_SRC})
TARGET_LINK_LIBRARIES(diffmtx cvtr kit)

ADD_EXECUTABLE(treedist ${TREEDIST_SRC})
TARGET_LINK_LIBRARIES(treedist cvtr kit)

ADD_EXECUTABLE(treeboot ${TREEBOOT_SRC})
TARGET_LINK_LIBRARIES(treeboot cvtr kit)

ADD_EXECUTABLE(dm2tree ${NJ_SRC})
TARGET_LINK_LIBRARIES(dm2tree cvtr kit)

ADD_EXECUTABLE(cvtree ${CVTREE_SRC})
TARGET_LINK_LIBRARIES(cvtree cvtr kit)

INSTALL(TARGETS g2cv cvdump cv2dm getdist 
  diffmtx treedist dm2tree cvtree mconv
  RUNTIME DESTINATION bin)
