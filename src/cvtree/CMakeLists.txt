# ~~~
# Copyright (c) 2026
# See the accompanying Manual for the contributors and the way to
# cite this work. Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2025-12-06 Saturday 14:35:00
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2026-05-26 Tuesday 14:46:03
# ~~~

include_directories(${CMAKE_SOURCE_DIR}/kit)
include_directories(${CMAKE_SOURCE_DIR}/cv)

SET(DMTRHEADS 
  distmatrix.h distmeth.h 
  tree.h treemeth.h 
  sampling.h sampleMeth.h
  )

SET(LIBDMTR_SRC ${DMTRHEADS}
  distmeth.cpp distmatrix.cpp 
  tree.cpp treemeth.cpp 
  sampling.cpp sampleMeth.cpp
  )

SET(CV_SRC ${DMTRHEADS} g2cv.h g2cv.cpp)

SET(CVDUMP_SRC ${DMTRHEADS} cvdump.cpp)

SET(DIST_SRC ${DMTRHEADS} cv2dm.h cv2dm.cpp)

SET(GETDIST_SRC ${DMTRHEADS} getdist.cpp)

SET(MCONV_SRC ${DMTRHEADS} mconv.cpp)

SET(DIFFMTX_SRC ${DMTRHEADS} diffMatrix.cpp)

SET(TREEDIST_SRC ${DMTRHEADS}
    treedist.h treedist.cpp
    setsym.h marktree.cpp marktree.h
  )

SET(TREEBOOT_SRC ${DMTRHEADS}
  treeboot.h treeboot.cpp
  setsym.h marktree.cpp marktree.h)

SET(NJ_SRC ${DMTRHEADS} dm2tree.h dm2tree.cpp)

SET(CVTREE_SRC ${DMTRHEADS}
  marktree.h marktree.cpp
  setsym.h cvtree.h cvtree.cpp)

ADD_LIBRARY(dmtr STATIC ${LIBDMTR_SRC})

ADD_EXECUTABLE(g2cv ${CV_SRC})
TARGET_LINK_LIBRARIES(g2cv dmtr cv kit)

ADD_EXECUTABLE(cvdump ${CVDUMP_SRC})
TARGET_LINK_LIBRARIES(cvdump dmtr cv kit)

ADD_EXECUTABLE(cv2dm ${DIST_SRC})
TARGET_LINK_LIBRARIES(cv2dm dmtr cv kit)

ADD_EXECUTABLE(getdist ${GETDIST_SRC})
TARGET_LINK_LIBRARIES(getdist dmtr cv kit)

ADD_EXECUTABLE(mconv ${MCONV_SRC})
TARGET_LINK_LIBRARIES(mconv dmtr cv kit)

ADD_EXECUTABLE(diffmtx ${DIFFMTX_SRC})
TARGET_LINK_LIBRARIES(diffmtx dmtr cv kit)

ADD_EXECUTABLE(treedist ${TREEDIST_SRC})
TARGET_LINK_LIBRARIES(treedist dmtr cv kit)

ADD_EXECUTABLE(treeboot ${TREEBOOT_SRC})
TARGET_LINK_LIBRARIES(treeboot dmtr cv kit)

ADD_EXECUTABLE(dm2tree ${NJ_SRC})
TARGET_LINK_LIBRARIES(dm2tree dmtr cv kit)

ADD_EXECUTABLE(cvtree ${CVTREE_SRC})
TARGET_LINK_LIBRARIES(cvtree dmtr cv kit)

INSTALL(TARGETS g2cv cvdump cv2dm getdist 
  diffmtx treedist dm2tree cvtree mconv
  RUNTIME DESTINATION bin)
