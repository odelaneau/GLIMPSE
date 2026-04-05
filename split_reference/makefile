#COMPILER MODE C++17
CXX?=g++


#create folders
dummy_build_folder_bin := $(shell mkdir -p bin)
dummy_build_folder_obj := $(shell mkdir -p obj)

#COMPILER & LINKER FLAGS
CXXFLAG=-std=c++17 -O3 -Wno-ignored-attributes -I../third_party/simde
LDFLAG=-O3

#CXXFLAG=-O0 -g -Wno-ignored-attributes
#LDFLAG=-O0

#DYNAMIC LIBRARIES
DYN_LIBS=-lz -lpthread -lbz2 -llzma -lcurl -lcrypto -ldeflate

HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

NAME=$(shell basename $(CURDIR))
BFILE=bin/GLIMPSE2_$(NAME)
EXEFILE=bin/GLIMPSE2_$(NAME)_static

#COMMIT_VERS=$(shell git rev-parse --short HEAD)
#COMMIT_DATE=$(shell git log -1 --format=%cd --date=short)
#CXXFLAG+= -D__COMMIT_ID__=\"$(COMMIT_VERS)\"
#CXXFLAG+= -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"

ARCH := $(shell uname -m)
ifeq ($(NAME),phase)
  ifneq (,$(filter x86_64 amd64,$(ARCH)))
    CXXFLAG+=-mavx2 -mfma
  endif
endif

COMMIT_VERS=3bed6d9
COMMIT_DATE=2022-12-07
ifneq ($(MAKECMDGOALS),clean)
 ifneq (, $(shell which git))
  COMMIT_VERS_GIT=$(shell git rev-parse --short HEAD 2>/dev/null)
  COMMIT_DATE_GIT=$(shell git log -1 --format=%cd --date=short 2>/dev/null)
  ifneq (, $(COMMIT_VERS_GIT))
   COMMIT_VERS=$(COMMIT_VERS_GIT)
  else
   $(info Cannot set commit version with git. Using default (last release) $(COMMIT_VERS))
  endif
  ifneq (, $(COMMIT_DATE_GIT))
   COMMIT_DATE=$(COMMIT_DATE_GIT)
  else
   $(info Cannot set commit date with git. Using default (last release) $(COMMIT_DATE))
  endif
 else
  $(info Git not available, cannot set commit date with git. Using default (last release)  $(COMMIT_VERS)/$(COMMIT_DATE))  
 endif
 $(info Commit version $(COMMIT_VERS))
 CXXFLAG+= -D__COMMIT_ID__=\"$(COMMIT_VERS)\" 
 $(info Commit date $(COMMIT_DATE))
 CXXFLAG+= -D__COMMIT_DATE__=\"$(COMMIT_DATE)\"
endif


# xSqueezeIt Support [YES/NO]
ifeq ($(XSI_SUPPORT),)
	XSI_SUPPORT=NO
endif
ifeq ($(XSI_SUPPORT),YES)
	XSI_INC=xsqueezeit/include
	XSI_LIB=xsqueezeit/libxsqueezeit.a
	DYN_LIBS+= -lzstd -lhts
	CXXFLAG+= -D__XSI__ -I$(XSI_INC)
endif

# BGEN Support [YES/NO]
ifeq ($(BGEN_SUPPORT),)
	BGEN_SUPPORT=NO
endif

ifeq ($(BGEN_SUPPORT),)
	BGEN_SUPPORT=NO
	BGEN_LIB=""
endif
ifeq ($(BGEN_SUPPORT),YES)
	BGEN_INC=../3rd_party/bgen/genfile/include/
	BGEN_LIB=../3rd_party/bgen/build/libbgen.a ../3rd_party/zstd/lib/libzstd.a
	#DYN_LIBS+=../3rd_party/zstd/lib/libzstd.a
	CXXFLAG+= -D__BGEN__ -I$(BGEN_INC)
endif

#CONDITIONAL PATH DEFINITON

# Auto-detected system target: works on linux/x86_64, linux/aarch64, macOS/arm64
OS := $(shell uname -s)
ifeq ($(OS),Darwin)
  SYS_CXX=clang++
  SYS_HTSSRC=/opt/homebrew
  SYS_HTSLIB_INC=$(SYS_HTSSRC)/include
  SYS_HTSLIB_LIB=$(SYS_HTSSRC)/lib/libhts.a
  SYS_BOOST_INC=/opt/homebrew/include
  SYS_BOOST_LIB=/opt/homebrew/lib
  SYS_DYN_LIBS=-L/opt/homebrew/lib -lz -lpthread -lbz2 -llzma -lcurl -lcrypto -ldeflate
else
  SYS_CXX=g++
  SYS_MULTIARCH=$(shell dpkg-architecture -qDEB_HOST_MULTIARCH 2>/dev/null || echo $(ARCH)-linux-gnu)
  SYS_HTSSRC=/usr/local
  SYS_HTSLIB_INC=$(SYS_HTSSRC)/include
  SYS_HTSLIB_LIB=$(SYS_HTSSRC)/lib/libhts.a
  SYS_BOOST_INC=/usr/include
  SYS_BOOST_LIB=/usr/lib/$(SYS_MULTIARCH)
  SYS_DYN_LIBS=-lz -lpthread -lbz2 -llzma -lcurl -lcrypto -ldeflate
endif

system: CXX=$(SYS_CXX)
system: HTSSRC=$(SYS_HTSSRC)
system: HTSLIB_INC=$(SYS_HTSLIB_INC)
system: HTSLIB_LIB=$(SYS_HTSLIB_LIB)
system: BOOST_INC=$(SYS_BOOST_INC)
system: BOOST_LIB_IO=$(SYS_BOOST_LIB)/libboost_iostreams.a
system: BOOST_LIB_PO=$(SYS_BOOST_LIB)/libboost_program_options.a
system: BOOST_LIB_SE=$(SYS_BOOST_LIB)/libboost_serialization.a
system: DYN_LIBS=$(SYS_DYN_LIBS)
system: $(BFILE)

desktop: HTSSRC=/home/srubinac/git
desktop: HTSLIB_INC=$(HTSSRC)/htslib
desktop: HTSLIB_LIB=$(HTSSRC)/htslib/libhts.a
desktop: BOOST_INC=/usr/include
desktop: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
desktop: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
desktop: BOOST_LIB_SE=/usr/lib/x86_64-linux-gnu/libboost_serialization.a
desktop: $(BFILE)

olivier: HTSSRC=$(HOME)/Tools
olivier: HTSLIB_INC=$(HTSSRC)/htslib-1.15
olivier: HTSLIB_LIB=$(HTSSRC)/htslib-1.15/libhts.a
olivier: BOOST_INC=/usr/include
olivier: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
olivier: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
olivier: BOOST_LIB_SE=/usr/lib/x86_64-linux-gnu/libboost_serialization.a
olivier: $(BFILE)

docker: HTSSRC=/usr/local
docker: HTSLIB_INC=$(HTSSRC)/include/htslib
docker: HTSLIB_LIB=$(HTSSRC)/lib/libhts.a
docker: BOOST_INC=/usr/include
docker: BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
docker: BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
docker: BOOST_LIB_SE=/usr/local/lib/libboost_serialization.a
docker: $(BFILE)

unil-dcsr: DYN_LIBS=-lz -lpthread -lcrypto $(XZ_ROOT)/lib/liblzma.so $(BZIP2_ROOT)/lib/libbz2.so $(CURL_ROOT)/lib/libcurl.so /dcsrsoft/spack/arolle/v1.0/spack/opt/spack/linux-rhel8-zen2/gcc-10.4.0/libdeflate-1.10-c5y7elzscnxjeh75hignyp23u5iprhl3/lib/libdeflate.a
unil-dcsr: HTSSRC=$(HTSLIB_ROOT)
unil-dcsr: HTSLIB_INC=$(HTSSRC)/include
unil-dcsr: HTSLIB_LIB=$(HTSSRC)/lib/libhts.a
unil-dcsr: BOOST_INC=$(BOOST_ROOT)/include
unil-dcsr: BOOST_LIB_IO=$(BOOST_ROOT)/lib/libboost_iostreams.a
unil-dcsr: BOOST_LIB_PO=$(BOOST_ROOT)/lib/libboost_program_options.a
unil-dcsr: BOOST_LIB_SE=$(BOOST_ROOT)/lib/libboost_serialization.a
unil-dcsr: $(BFILE)

static_exe: HTSSRC=/home/srubinac/git
static_exe: HTSLIB_INC=$(HTSSRC)/htslib_minimal
static_exe: HTSLIB_LIB=$(HTSSRC)/htslib_minimal/libhts.a
static_exe: BOOST_INC=/usr/include
static_exe: BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
static_exe: BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
static_exe: BOOST_LIB_SE=/usr/local/lib/libboost_serialization.a

static_exe: $(EXEFILE)

#COMPILATION RULES
all: desktop
	
$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) $(BOOST_LIB_SE) $(BGEN_LIB) -o $@ $(DYN_LIBS)

$(EXEFILE): $(OFILE)
	$(CXX) $(LDFLAG) -static -static-libgcc -static-libstdc++ -pthread -o $(EXEFILE) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) $(BOOST_LIB_SE) -Wl,-Bstatic $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE)
