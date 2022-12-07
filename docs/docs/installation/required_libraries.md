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
GLIMPSE2 requires several libraries installed on the system. Here we assume most of the libraries are available, and we focus on two main libraries:

- HTSlib version >= 1.7: A C library for reading/writing high-throughput sequencing data.
- BOOST version >= 1.65: A set of peer-reviewed portable C++ source libraries. GLIMPSE2 uses two specific BOOST libraries: `iostreams` and `program_options`.

### HTSlib
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


### Boost
As GLIMPSE2 only requires few of the boost libraries, we'll build the smallest possible boost build, without requiring root privileges. Please refer to the [Boost installation instructions](https://www.boost.org/doc/libs/1_73_0/more/getting_started/unix-variants.html#easy-build-and-install) for complete details. Here we provide a basic script to the minimal build of Boost v1.73.0 required to run GLIMPSE2:

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

After this, you will also find a copy of the Boost headers in the include/ subdirectory of the installation prefix (our current directory). The Boost static libraries (`libboost_iostreams.a` and `libboost_program_options.a`) are found in the subfolder `./lib/` of your installation prefix.

### Additional libraries

Make sure that the following standard library flags can be used by g++ on your system:
- `-lz`,`-lbz2` and `-llzma` for reading/writing compressed files.
- `-lm` for basic math operations.
- `-lpthread` for multi-threading

You can do so by checking the outcome of the following commands:
<div class="code-example" markdown="1">
```bash
locate -b '\libz.so'
locate -b '\libbz2.so'
locate -b '\liblzma.so'
locate -b '\libm.so'
locate -b '\libpthread.so'
locate -b '\libcurl.so'
```
</div>


