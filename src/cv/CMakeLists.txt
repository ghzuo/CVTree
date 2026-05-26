# ~~~
# Copyright (c) 2026
# See the accompanying Manual for the contributors and the way to
# cite this work. Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2025-12-06 Saturday 14:35:00
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2026-05-26 Tuesday 14:47:11
# ~~~

include_directories(${CMAKE_SOURCE_DIR}/kit)

SET(CVHEADS 
  genome.h kstring.h 
  cv.h cvbundle.h aligncv.h karray.h cvarray.h cvmeth.h
  triangleMatrix.h rectangleMatrix.h similarity.h
)

SET(LIBCV_SRC ${CVTREEHEADS}
  genome.cpp kstring.cpp 
  cv.cpp cvbundle.cpp aligncv.cpp karray.cpp cvarray.cpp cvmeth.cpp
  triangleMatrix.cpp rectangleMatrix.cpp similarity.cpp
)

ADD_LIBRARY(cv STATIC ${LIBCV_SRC})

