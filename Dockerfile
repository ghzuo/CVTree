## Stage for build cvtree
FROM ubuntu AS dev
LABEL Version=0.1 \
  MAINTAINER="Guanghong Zuo<ghzuo@fudan.edu.cn>"\
  description="Docker image for CVTree" 

## for develop environment
RUN apt-get update -yqq
RUN apt-get install -yqq build-essential
RUN apt-get install -yqq cmake
RUN apt-get install -yqq libhdf5-dev h5utils

## Build cvtree
WORKDIR /root
COPY ./src /root/cvtree/src
COPY ./CMakeLists.txt /root/cvtree/
RUN mkdir cvtree/build/ && cd cvtree/build/ && cmake .. && make 

## Stage for run cvtree 
FROM ubuntu AS run
COPY --from=dev /root/cvtree/build/bin/* /usr/local/bin/
RUN apt-get update -yqq && apt-get install -yqq libgomp1 libsz2 && rm -rf /var/lib/apt/lists/* && apt-get clean

## for workplace
WORKDIR /root
CMD [ "/bin/bash" ]