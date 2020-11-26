###
# Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai, China.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2017-09-01 15:55:03
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2020-11-26 09:17:46
###

LINK_LIBRARIES(libkit.a)

SET(KITHEADS 
  stringOpt.h  memory.h
  info.h
  )

SET(LIBKIT_SRC ${KITHEADS}
  stringOpt.cpp  memory.cpp
  info.cpp 
  )

ADD_LIBRARY(kit STATIC ${LIBKIT_SRC})

