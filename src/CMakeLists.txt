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

SET(NJ_SRC ${CVTREEHEADS}
  neighbor.h neighbor.cpp)

SET(CVTREE_SRC ${CVTREEHEADS}
  cvtree.h cvtree.cpp)

ADD_LIBRARY(cvtr STATIC ${LIBCVTREE_SRC})

ADD_EXECUTABLE(cv ${CV_SRC})
TARGET_LINK_LIBRARIES(cv cvtr)

ADD_EXECUTABLE(cvdump ${CVDUMP_SRC})
TARGET_LINK_LIBRARIES(cvdump cvtr)

ADD_EXECUTABLE(dist ${DIST_SRC})
TARGET_LINK_LIBRARIES(dist cvtr)

ADD_EXECUTABLE(getdist ${GETDIST_SRC})
TARGET_LINK_LIBRARIES(getdist cvtr)

ADD_EXECUTABLE(nj ${NJ_SRC})
TARGET_LINK_LIBRARIES(nj cvtr)

ADD_EXECUTABLE(cvtree ${CVTREE_SRC})
TARGET_LINK_LIBRARIES(cvtree cvtr)

INSTALL(TARGETS cv cvdump dist getdist nj cvtree
   RUNTIME DESTINATION bin)
