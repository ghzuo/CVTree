###
# Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai, China.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2017-09-01 15:55:03
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2018-04-26 19:51:32
###

LINK_LIBRARIES(libcvtr.a)

SET(CVTREEHEADS 
  readgenome.h kstring.h
  stringOpt.h  memory.h
  distmatrix.h info.h
  method.h neighborJoint.h
  distance.h)

SET(LIBCVTREE_SRC ${CVTREEHEADS}
  readgenome.cpp kstring.cpp
  stringOpt.cpp  memory.cpp
  distmatrix.cpp tree.cpp
  info.cpp method.cpp
  neighborJoint.cpp distance.cpp
  )

SET(CV_SRC ${CVTREEHEADS}
  cv.h cv.cpp)

SET(CVDUMP_SRC ${CVTREEHEADS}
  cvdump.cpp)

SET(DIST_SRC ${CVTREEHEADS}
  dist.h dist.cpp)

SET(GETDIST_SRC ${CVTREEHEADS}
  getdist.cpp)

SET(DIFFMTX_SRC ${CVTREEHEADS} diffMatrix.cpp)

SET(NJ_SRC ${CVTREEHEADS}
  neighbor.h neighbor.cpp)

SET(CVTREE_SRC ${CVTREEHEADS}
  cvtree.h cvtree.cpp)

ADD_LIBRARY(cvtr STATIC ${LIBCVTREE_SRC})

ADD_EXECUTABLE(g2cv ${CV_SRC})
TARGET_LINK_LIBRARIES(g2cv cvtr)

ADD_EXECUTABLE(cvdump ${CVDUMP_SRC})
TARGET_LINK_LIBRARIES(cvdump cvtr)

ADD_EXECUTABLE(cv2dm ${DIST_SRC})
TARGET_LINK_LIBRARIES(cv2dm cvtr)

ADD_EXECUTABLE(getdist ${GETDIST_SRC})
TARGET_LINK_LIBRARIES(getdist cvtr)

ADD_EXECUTABLE(diffmtx ${DIFFMTX_SRC})
TARGET_LINK_LIBRARIES(diffmtx cvtr)

ADD_EXECUTABLE(dm2tree ${NJ_SRC})
TARGET_LINK_LIBRARIES(dm2tree cvtr)

ADD_EXECUTABLE(cvtree ${CVTREE_SRC})
TARGET_LINK_LIBRARIES(cvtree cvtr)

INSTALL(TARGETS g2cv cvdump cv2dm getdist diffmtx dm2tree cvtree
   RUNTIME DESTINATION bin)
