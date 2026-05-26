# ~~~
# Copyright (c) 2026
# See the accompanying Manual for the contributors and the way to
# cite this work. Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2025-12-06 Saturday 14:35:00
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2026-05-26 Tuesday 14:48:40
# ~~~

SET(KITHEADS 
  stringOpt.h  memory.h jsonctl.h
  info.h ompOpt.h kit.h
  )

SET(LIBKIT_SRC ${KITHEADS}
  stringOpt.cpp  memory.cpp
  info.cpp fileOpt.cpp jsonctl.cpp ompOpt.cpp
  )

ADD_LIBRARY(kit STATIC ${LIBKIT_SRC})
