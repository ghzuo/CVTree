###
# Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai, China.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2020-10-12 14:35:55
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: Fri Jul 12 2024
###

## Stage for build cvtree
FROM alpine AS dev
LABEL Version=0.2 \
  MAINTAINER="Guanghong Zuo<ghzuo@ucas.ac.cn>"\
  description="Docker image for CVTree" 

## for develop environment
RUN apk --update add --no-cache g++ make cmake zlib-dev boost-dev
RUN apk --update add --no-cache hdf5-dev --repository=http://dl-cdn.alpinelinux.org/alpine/edge/community
RUN apk --update add --no-cache nlohmann-json --repository=http://dl-cdn.alpinelinux.org/alpine/edge/community


## Build cvtree
WORKDIR /root
COPY ./cvtree /root/cvtree/cvtree
COPY ./kit /root/cvtree/kit
COPY ./CMakeLists.txt /root/cvtree/
RUN mkdir cvtree/build/ && cd cvtree/build/ && cmake .. && make 

## Stage for run cvtree 
FROM alpine AS run
COPY --from=dev /root/cvtree/build/bin/* /usr/local/bin/
RUN apk --update add --no-cache libgomp libstdc++
RUN apk --update add --no-cache hdf5-dev --repository=http://dl-cdn.alpinelinux.org/alpine/edge/community

## for workplace
WORKDIR /root/data
ENTRYPOINT ["cvtree"]
