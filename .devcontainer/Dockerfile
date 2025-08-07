FROM ubuntu:20.04

LABEL org.opencontainers.image.created="2022-11-30"
LABEL org.opencontainers.image.url="https://github.com/odelaneau/GLIMPSE"
LABEL org.opencontainers.image.version="2.0.0"
LABEL org.opencontainers.image.licences="MIT"
LABEL org.opencontainers.image.title="glimpse"
LABEL org.opencontainers.image.authors="simone.rubinacci@unil.ch"

WORKDIR /docker_build/
# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York

# Install required packages for build
RUN apt-get update && apt-get install -y build-essential=12.8ubuntu1.1 \
                                         libbz2-dev=1.0.8-2 \
                                         libcurl4-openssl-dev=7.68.0-1ubuntu2.25 \
                                         autoconf=2.69-11.1 \
                                         libssl-dev=1.1.1f-1ubuntu2.24 \
                                         wget=1.20.3-1ubuntu2.1 \
                                         zlib1g-dev=1:1.2.11.dfsg-2ubuntu1.5 \
                                         liblzma-dev=5.2.4-1ubuntu1.1 \
                                         libdeflate-dev=1.5-3 \
                                         python3=3.8.2-0ubuntu2 \
                                         libboost-regex-dev=1.71.0.0ubuntu2\
                                         libboost-iostreams-dev=1.71.0.0ubuntu2\
                                         libboost-program-options-dev=1.71.0.0ubuntu2\
                                         libhts-dev=1.10.2-3ubuntu0.1\
                                         libboost-serialization-dev=1.71.0.0ubuntu2

# Install required packages for docs 
RUN apt-get update && apt-get install -y ruby-full=1:2.7+1
RUN gem install bundler -v 2.3.27

# Download and build boost program_options and iostreams
RUN wget https://archives.boost.io/release/1.78.0/source/boost_1_78_0.tar.gz && \
tar -xf boost_1_78_0.tar.gz && \
rm boost_1_78_0.tar.gz && \
cd boost_1_78_0/ && \
./bootstrap.sh --with-libraries=iostreams,program_options,serialization --prefix=../boost && \
./b2 install && \
cd .. && \
cp boost/lib/libboost_iostreams.a boost/lib/libboost_program_options.a boost/lib/libboost_serialization.a /usr/local/lib/ && \
cp -r boost/include/boost/ /usr/include/ && \
rm -r boost_1_78_0 boost

# Download and build htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
tar -xf htslib-1.16.tar.bz2 && \
rm htslib-1.16.tar.bz2 && \
cd htslib-1.16 && \
autoheader && \
autoconf && \
./configure --enable-libcurl && \
make install && \
cd .. && \
rm -r htslib-1.16

# Install github
RUN apt-get install -y --no-install-recommends git=1:2.25.1-1ubuntu3.14

# Install cmake
RUN mkdir -p /usr/local
RUN cd /usr/local &&\
    wget https://github.com/Kitware/CMake/releases/download/v4.0.2/cmake-4.0.2-linux-x86_64.sh &&\
    chmod +x cmake-4.0.2-linux-x86_64.sh &&\
    ./cmake-4.0.2-linux-x86_64.sh --skip-license

# Install required packages to debug
ARG INSTALL_DEBUG_TOOLS=true
RUN if [ "$INSTALL_DEBUG_TOOLS" = "true" ]; then \
    apt-get install -y --no-install-recommends openssh-client=1:8.2p1-4ubuntu0.13 &&\
    apt-get install -y --no-install-recommends gdb=9.2-0ubuntu1~20.04.2 &&\
    apt-get install -y --no-install-recommends bcftools=1.10.2-2 &&\
    apt-get install -y --no-install-recommends less=551-1ubuntu0.3 &&\
    apt-get install -y --no-install-recommends samtools=1.10-3 &&\
    apt-get install -y --no-install-recommends bash-completion=1:2.10-1ubuntu1 ;\
    fi

# Set up terminal 
COPY .devcontainer/.bashrc /root/.bashrc


# Set workspace
RUN mkdir -p /workspace
WORKDIR /workspace
