## Stage for build cvtree
FROM alpine AS dev
LABEL Version=0.1 \
  MAINTAINER="Guanghong Zuo<ghzuo@fudan.edu.cn>"\
  description="Docker image for CVTree" 

## for develop environment
RUN apk --update add --no-cache g++ make cmake zlib-dev
RUN apk --update add --no-cache hdf5-dev hdf5-static --repository=http://dl-cdn.alpinelinux.org/alpine/edge/community


## Build cvtree
WORKDIR /root
COPY ./src /root/cvtree/src
COPY ./CMakeLists.txt /root/cvtree/
RUN mkdir cvtree/build/ && cd cvtree/build/ && cmake .. && make 

## Stage for run cvtree 
FROM alpine AS run
COPY --from=dev /root/cvtree/build/bin/* /usr/local/bin/
RUN apk --update add --no-cache libgomp libstdc++

## for workplace
WORKDIR /root/data
