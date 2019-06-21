FROM ubuntu
LABEL Version=0.1 \
  MAINTAINER="Guanghong Zuo<ghzuo@fudan.edu.cn>"\
  description="Docker image for CVTree" 

## for develop environment
RUN apt-get update -yqq
RUN apt-get install -yqq build-essential
RUN apt-get install -yqq cmake
RUN apt-get install -yqq libhdf5-dev h5utils
RUN apt-get install -yqq libnetcdf-cxx-legacy-dev netcdf-bin nco

## Build cvtree
WORKDIR /root
COPY ./src /root/cvtree/src
COPY ./CMakeLists.txt /root/cvtree/
RUN mkdir cvtree/build/ && cd cvtree/build/ && cmake .. && make install 
RUN rm -rf /root/cvtree

## clean apt-get
RUN apt-get clean

## for workplace
CMD [ "/bin/bash" ]