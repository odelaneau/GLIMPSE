---
layout: default
title: Required libraries
nav_order: 2
parent: Build from source
grand_parent: Installation
permalink: /docs/installation/build_from_source/required_libraries
---
# Required libraries
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Required libraries
GLIMPSE2 requires several libraries installed on the system:

- **HTSlib** version >= 1.7: A C library for reading/writing high-throughput sequencing data.
- **Boost** version >= 1.65: A set of peer-reviewed portable C++ source libraries. GLIMPSE2 uses three specific Boost libraries: `iostreams`, `program_options` and `serialization`.
- **SIMDe** (bundled): The [SIMDe](https://github.com/simd-everywhere/simde) header-only library is included as a git submodule and provides portable SIMD support across x86_64 and ARM64 platforms.

### macOS (Homebrew)

On macOS with [Homebrew](https://brew.sh/), all required libraries can be installed with a single command:

<div class="code-example" markdown="1">
```bash
brew install htslib boost libdeflate
```
</div>

### Linux (package manager)

Most Linux distributions package all required libraries. Installing from your distribution's package manager is the easiest approach.

#### Debian / Ubuntu
<div class="code-example" markdown="1">
```bash
sudo apt-get install build-essential pkg-config \
  libhts-dev libhtscodecs-dev libboost-iostreams-dev libboost-program-options-dev \
  libboost-serialization-dev libdeflate-dev zlib1g-dev libbz2-dev \
  liblzma-dev libssl-dev
```
</div>

#### Fedora / RHEL / Rocky Linux
<div class="code-example" markdown="1">
```bash
sudo dnf install gcc-c++ make pkg-config \
  htslib-devel boost-devel boost-static libdeflate-devel zlib-devel bzip2-devel \
  xz-devel libcurl-devel openssl-devel
```
</div>

Not all Fedora-based distributions (e.g., Amazon Linux) make `htslib-devel` and `libdeflate-devel` available in their package repositories. If these packages are not found, you will need to build them from source (see below).

After installing packages, build with `make system` from the GLIMPSE2 root directory. The build system will automatically locate libraries using `pkg-config` (for HTSlib) and by searching standard installation paths.

### Linux (from source)

If your distribution does not package these libraries, or you need specific versions, you can build them from source.

#### HTSlib
Building HTSlib is straightforward and does not require root privileges. Please refer to the [HTSlib](http://www.htslib.org/) documentation for complete details. Here we provide a basic script to install HTSlib v1.16:

<div class="code-example" markdown="1">
```bash
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xf htslib-1.16.tar.bz2
mv htslib-1.16 htslib
cd htslib
autoheader; autoconf; ./configure; #optional
make
```
</div>

After this, you'll find the libhts.a inside your current directory and the include files inside subdirectory: `./include/`


#### Boost
As GLIMPSE2 only requires few of the Boost libraries, we'll build the smallest possible Boost build, without requiring root privileges. Please refer to the [Boost installation instructions](https://www.boost.org/doc/libs/1_73_0/more/getting_started/unix-variants.html#easy-build-and-install) for complete details. Here we provide a basic script to the minimal build of Boost v1.73.0 required to run GLIMPSE2:

<div class="code-example" markdown="1">
```bash
wget https://boostorg.jfrog.io/artifactory/main/release/1.73.0/source/boost_1_73_0.tar.bz2
tar --bzip2 -xf boost_1_73_0.tar.bz2
cd boost_1_73_0
./bootstrap.sh --with-libraries=iostreams,program_options,serialization --prefix=../boost #where ../boost is your custom boost installation prefix
./b2 install
cd ../boost #change this to the folder you used as --prefix for the bootstrap script
```
</div>

After this, you will also find a copy of the Boost headers in the include/ subdirectory of the installation prefix (our current directory). The Boost static libraries (`libboost_iostreams.a`, `libboost_program_options.a` and `libboost_serialization.a`) are found in the subfolder `./lib/` of your installation prefix.

#### Additional libraries

Make sure that the following standard library flags can be used by your compiler on your system:
- `-lz`,`-lbz2` and `-llzma` for reading/writing compressed files.
- `-lm` for basic math operations.
- `-lpthread` for multi-threading
- `-lcurl` for network access.
- `-lcrypto` for cryptographic functions.
- `-ldeflate` for libdeflate compression.


