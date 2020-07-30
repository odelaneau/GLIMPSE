cmake_minimum_required(VERSION 3.5) 

project(GLIMPSE
  VERSION 1.0.0
  LANGUAGES CXX
  DESCRIPTION "Low Coverage Calling of Genotypes"
  HOMEPAGE_URL "https://odelaneau.github.io/GLIMPSE/")

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo"
    FORCE)
endif ()

option(STATIC_BUILD "Build static executation" OFF)

if (STATIC_BUILD)
  set(CMAKE_FIND_LIBRARY_PREFIXES ".a")
endif (STATIC_BUILD)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Static build: ${STATIC_BUILD}")

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost COMPONENTS iostreams program_options)
find_package(HTSlib REQUIRED)
find_package(CURL REQUIRED)
find_package(OpenSSL REQUIRED)

if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif ()

set(CMAKE_CXX_FLAGS_DEBUG "-g --ggdb -O1 -D_LIBCPP_DEBUG -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

if (STATIC_BUILD)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl, --whole-archive -lpthread -Wl,--no-whole-archive -static -static-libgcc -static-libstdc++")
endif (STATIC_BUILD)

add_subdirectory(chunk)
add_subdirectory(concordance)
add_subdirectory(ligate)
add_subdirectory(phase)
#add_subdirectory(rephase)
add_subdirectory(sample)
add_subdirectory(snparray)
add_subdirectory(stats)
