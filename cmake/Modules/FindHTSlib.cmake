# - Try to find HTSlib
# Once done, this will define
#
#  HTSlib_FOUND - system has HTSlib
#  HTSlib_INCLUDE_DIRS - the HTSlib include directories
#  HTSlib_LIBRARIES - link these to use HTSlib

include(LibFindMacros)

# HTSlib dependencies
libfind_package(HTSLib ZLIB REQUIRED)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(HTSlib_PKGCONF HTSlib)

# Locate include dir
find_path(HTSlib_INCLUDE_DIR
  NAMES hts.h sam.h
  PATHS ${HTSlib_PKGCONF_INCLUDE_DIRS}
  PATH_SUFFIXES include include/htslib
)

# Find the library itself
find_library(HTSlib_LIBRARY
  NAMES hts
  PATHS ${HTSlib_PKGCONF_LIBRARY_DIRS}
)

libfind_process(HTSlib)
