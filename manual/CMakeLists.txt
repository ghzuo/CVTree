###
# Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai, China.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2018-08-13 21:29:23
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2018-08-13 21:29:23
###

include(UseLATEX.cmake)

ADD_LATEX_DOCUMENT(
  cvtreeManual.tex
  # BIBFILES bibliography.bib
  IMAGE_DIRS fig
  TARGET_NAME manual
  EXCLUDE_FROM_ALL
)

# explicitly add dependency to doc
# add_dependencies(doc manual)

add_custom_command(TARGET manual POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_BINARY_DIR}/cvtreeManual.pdf
      ${CMAKE_BINARY_DIR}/cvtreeMannal.pdf)