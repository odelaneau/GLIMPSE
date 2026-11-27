---
title: "Build the dependencies yourself"
description: "Build HTSlib and Boost in your home directory, without root access."
---

[Building GLIMPSE2](../build/) normally uses the libraries from your package manager. Build HTSlib and Boost yourself if you have no root access, if your distribution doesn't package them or ships versions that are too old, or if you want Boost 1.84 or newer so that you can read and write `.bin` reference panels that other builds can use (see the note on the [build page](../build/)).

The commands below install both libraries into one directory, `$DEPS`, that you own. Set it to any location you like.

```bash
export DEPS=$HOME/glimpse-deps
mkdir -p $DEPS/src && cd $DEPS/src

# HTSlib 1.25
wget https://github.com/samtools/htslib/releases/download/1.25/htslib-1.25.tar.bz2
tar -xf htslib-1.25.tar.bz2
cd htslib-1.25
./configure --prefix=$DEPS
make -j4 && make install
cd ..

# Boost 1.90, only the three libraries GLIMPSE2 uses
wget https://archives.boost.io/release/1.90.0/source/boost_1_90_0.tar.bz2
tar --bzip2 -xf boost_1_90_0.tar.bz2
cd boost_1_90_0
./bootstrap.sh --with-libraries=iostreams,program_options,serialization --prefix=$DEPS
./b2 install
```

HTSlib needs zlib, bzip2, xz, libcurl, OpenSSL and libdeflate to be installed on the system. Its [documentation](http://www.htslib.org/) lists the details. Compiling Boost takes a few minutes.

After this, `$DEPS/include` holds the headers of both libraries and `$DEPS/lib` their static libraries. Build GLIMPSE2 against them:

```bash
cd /path/to/GLIMPSE
make system SYS_HTSLIB_INC=$DEPS/include SYS_HTSLIB_LIBDIR=$DEPS/lib \
  SYS_BOOST_INC=$DEPS/include SYS_BOOST_LIB=$DEPS/lib
```
