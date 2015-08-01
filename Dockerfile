FROM ubuntu:14.04

# Add openMVG binaries to path
ENV PATH $PATH:/opt/openMVG_Build/software/SfM

# Get dependencies
RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  graphviz \
  git \
  gcc-4.6 \ 
  g++-4.6 \ 
  gcc-4.6-multilib \  
  g++-4.6-multilib \ 
  libpng-dev \
  libjpeg-dev \
  libtiff-dev \
  libxxf86vm1 \
  libxxf86vm-dev \
  libxi-dev \
  libxrandr-dev \
  python-dev \  
  python-pip

# Clone the openvMVG repo 
RUN cd /opt && git clone --recursive https://github.com/openMVG/openMVG.git

# Build
RUN mkdir /opt/openMVG_Build && cd /opt/openMVG_Build && cmake -DCMAKE_BUILD_TYPE=RELEASE \
  -DOpenMVG_BUILD_TESTS=ON -DOpenMVG_BUILD_EXAMPLES=ON . ../openMVG/src/ && make

RUN cd /opt/openMVG_Build && make test
