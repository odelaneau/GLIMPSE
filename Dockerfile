# syntax=docker/dockerfile:1
#
# GLIMPSE2 tools, bcftools, bgzip and tabix for linux/amd64 and linux/arm64.
# Build with ./docker-build.sh [ref] [-t tag] [--push repository], which sets the build args.
# By hand, from a checkout with the SIMDe submodule:
#
#   docker build -t glimpse:2.1.0 \
#     --build-arg GLIMPSE_COMMIT=$(git rev-parse --short HEAD) \
#     --build-arg GLIMPSE_DATE=$(git log -1 --format=%cd --date=short) .
#
# amd64 needs AVX2 and FMA. See docs/docs/installation/docker.md.

ARG UBUNTU_VERSION=26.04
ARG GCC_VERSION=15
ARG LLVM_VERSION=23

FROM ubuntu:${UBUNTU_VERSION} AS build-base
ARG GCC_VERSION
ARG LLVM_VERSION
# Boost 1.90 writes .bin archive version 20, readable by Boost 1.84 and later.
# clang builds much faster SIMDe NEON kernels than g++ on arm64.
# bsdtar, because 26.04's GNU tar cannot create directories under Rosetta.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential g++-${GCC_VERSION} ca-certificates wget bzip2 libarchive-tools pkgconf \
      zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libdeflate-dev \
      libboost-iostreams-dev libboost-program-options-dev libboost-serialization-dev \
 && if [ "$(dpkg --print-architecture)" = arm64 ]; then \
      . /etc/os-release \
      && wget -qO /etc/apt/keyrings/apt.llvm.org.asc https://apt.llvm.org/llvm-snapshot.gpg.key \
      && echo "deb [signed-by=/etc/apt/keyrings/apt.llvm.org.asc] https://apt.llvm.org/${VERSION_CODENAME}/ llvm-toolchain-${VERSION_CODENAME}-${LLVM_VERSION} main" \
           > /etc/apt/sources.list.d/apt.llvm.org.list \
      && apt-get update \
      && apt-get install -y --no-install-recommends clang-${LLVM_VERSION}; \
    fi \
 && rm -rf /var/lib/apt/lists/*

# Ubuntu 26.04 ships HTSlib 1.22; GLIMPSE2_concordance can crash with 1.22 and 1.23.
FROM build-base AS htslib
ARG HTSLIB_VERSION=1.25
ARG HTSLIB_SHA256=158c271e4e2986bb2dd11283b640edc10c56695b293df77a109ce99480303730
WORKDIR /build
RUN wget -q https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
 && echo "${HTSLIB_SHA256}  htslib-${HTSLIB_VERSION}.tar.bz2" | sha256sum -c - \
 && bsdtar -xf htslib-${HTSLIB_VERSION}.tar.bz2 \
 && cd htslib-${HTSLIB_VERSION} \
 && ./configure --prefix=/usr/local --enable-libcurl --enable-gcs --enable-s3 \
 && make -j"$(nproc)" \
 && make install \
 && mkdir -p /runtime/bin /runtime/lib \
 && cp -a /usr/local/lib/libhts.so* /runtime/lib/ \
 && cp /usr/local/bin/bgzip /usr/local/bin/tabix /usr/local/bin/htsfile /runtime/bin/ \
 && strip /runtime/bin/* /runtime/lib/libhts.so.*.*

FROM htslib AS bcftools
ARG BCFTOOLS_VERSION=1.25
ARG BCFTOOLS_SHA256=f96e583d2315998645bb56756f1c8f230c489512db2655b0f623198ebc8a2d99
RUN wget -q https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
 && echo "${BCFTOOLS_SHA256}  bcftools-${BCFTOOLS_VERSION}.tar.bz2" | sha256sum -c - \
 && bsdtar -xf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
 && cd bcftools-${BCFTOOLS_VERSION} \
 && ./configure --prefix=/usr/local --with-htslib=system \
 && make -j"$(nproc)" \
 && make install \
 && rm -rf /runtime && mkdir -p /runtime/bin /runtime/libexec \
 && cp /usr/local/bin/bcftools /runtime/bin/ \
 && cp -a /usr/local/libexec/bcftools /runtime/libexec/ \
 && strip /runtime/bin/bcftools /runtime/libexec/bcftools/*.so

FROM htslib AS glimpse
ARG GCC_VERSION
ARG LLVM_VERSION
# Commit and date the tools report (the build has no git history).
ARG GLIMPSE_COMMIT=unknown
ARG GLIMPSE_DATE=unknown
COPY . /src/GLIMPSE
WORKDIR /src/GLIMPSE
RUN test -f third_party/simde/simde/x86/avx2.h \
      || { echo "third_party/simde is empty: run 'git submodule update --init' and build again" >&2; exit 1; } \
 && if [ "$(dpkg --print-architecture)" = arm64 ]; then CXX=clang++-${LLVM_VERSION}; else CXX=g++-${GCC_VERSION}; fi \
 && make system -j"$(nproc)" CXX="$CXX" COMMIT_VERS="$GLIMPSE_COMMIT" COMMIT_DATE="$GLIMPSE_DATE" \
 && rm -rf /runtime && mkdir -p /runtime/bin \
 && cp */bin/GLIMPSE2_* /runtime/bin/ \
 && strip /runtime/bin/*

FROM ubuntu:${UBUNTU_VERSION}
ARG GLIMPSE_COMMIT=unknown
LABEL org.opencontainers.image.title="GLIMPSE2" \
      org.opencontainers.image.description="Low-coverage imputation and phasing (GLIMPSE2 tools, bcftools, bgzip, tabix)" \
      org.opencontainers.image.source="https://github.com/odelaneau/GLIMPSE" \
      org.opencontainers.image.url="https://odelaneau.github.io/GLIMPSE" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.revision="${GLIMPSE_COMMIT}"
# ca-certificates for https, S3 and GCS URLs.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      ca-certificates libcurl4t64 libssl3t64 libdeflate0 liblzma5 libbz2-1.0 zlib1g \
 && rm -rf /var/lib/apt/lists/*
# Copying the /runtime directories keeps symlinks.
COPY --from=htslib /runtime/ /usr/local/
COPY --from=bcftools /runtime/ /usr/local/
COPY --from=glimpse /runtime/ /usr/local/
RUN ldconfig \
 && if ldd /usr/local/bin/* /usr/local/libexec/bcftools/*.so | grep "not found"; then exit 1; fi \
 && for tool in chunk split_reference phase ligate concordance inspect; do GLIMPSE2_$tool --help > /dev/null; done \
 && bcftools --version | head -1 \
 && bgzip --version | head -1
