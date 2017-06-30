FROM ubuntu:14.04

# Add openMVG binaries to path
ENV PATH $PATH:/opt/openMVG_Build/install/bin

# Get dependencies
RUN apt-get update && apt-get install -y \
  build-essential \
  graphviz \
  git \
  gcc-4.8 \ 
  gcc-4.8-multilib \  
  libpng-dev \
  libjpeg-dev \
  libtiff-dev \
  libxxf86vm1 \
  libxxf86vm-dev \
  libxi-dev \
  libxrandr-dev \
  python-dev \  
  python-pip

# build cmake (ubuntu 14.04 comes with cmake 2.8, we want a 3.X)
RUN apt-get install -y curl
RUN curl -O https://cmake.org/files/v3.8/cmake-3.8.0.tar.gz \
     && tar -xvf cmake-3.8.0.tar.gz
RUN cd cmake-3.8.0 && ./bootstrap && make && make install

# Clone the openvMVG repo 
ADD . /opt/openMVG
RUN cd /opt/openMVG && git submodule update --init --recursive

# Build
RUN mkdir /opt/openMVG_Build && cd /opt/openMVG_Build && cmake -DCMAKE_BUILD_TYPE=RELEASE \
  -DCMAKE_INSTALL_PREFIX="/opt/openMVG_Build/install" -DOpenMVG_BUILD_TESTS=ON \
  -DOpenMVG_BUILD_EXAMPLES=ON . ../openMVG/src/ && make -j 4

RUN cd /opt/openMVG_Build && make test
