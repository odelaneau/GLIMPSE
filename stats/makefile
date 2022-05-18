#COMPILER MODE C++11
CXX=g++ -std=c++11


#create folders
dummy_build_folder_bin := $(shell mkdir -p bin)
dummy_build_folder_obj := $(shell mkdir -p obj)

#COMPILER & LINKER FLAGS
CXXFLAG=-O3
LDFLAG=-O3

#DYNAMIC LIBRARIES
DYN_LIBS=-lz -lpthread -lbz2 -llzma -lcurl -lcrypto

HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

NAME=$(shell basename $(CURDIR))
BFILE=bin/GLIMPSE_$(NAME)
EXEFILE=bin/GLIMPSE_$(NAME)_static

#CONDITIONAL PATH DEFINITON
system: DYN_LIBS=-lz -lpthread -lbz2 -llzma
system: HTSSRC=/usr/local
system: HTSLIB_INC=$(HTSSRC)/include/htslib
system: HTSLIB_LIB=$(HTSSRC)/lib/libhts.a
system: BOOST_INC=/usr/include
system: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
system: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
system: BOOST_LIB_SE=/usr/lib/x86_64-linux-gnu/libboost_serialization.a
system: $(BFILE)

desktop: HTSSRC=../..
desktop: HTSLIB_INC=$(HTSSRC)/htslib
desktop: HTSLIB_LIB=$(HTSSRC)/htslib/libhts.a
desktop: BOOST_INC=/usr/include
desktop: BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
desktop: BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
desktop: $(BFILE)

olivier: HTSSRC=$(HOME)/Tools
olivier: HTSLIB_INC=$(HTSSRC)/htslib-1.9
olivier: HTSLIB_LIB=$(HTSSRC)/htslib-1.9/libhts.a
olivier: BOOST_INC=/usr/include
olivier: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
olivier: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
olivier: $(BFILE)

unil-dcsr: DYN_LIBS=-lz -lpthread -lcrypto $(XZ_ROOT)/lib/liblzma.so $(BZIP2_ROOT)/lib/libbz2.so $(CURL_ROOT)/lib/libcurl.so
unil-dcsr: HTSSRC=$(HTSLIB_ROOT)
unil-dcsr: HTSLIB_INC=$(HTSSRC)/include
unil-dcsr: HTSLIB_LIB=$(HTSSRC)/lib/libhts.a
unil-dcsr: BOOST_INC=$(BOOST_ROOT)/include
unil-dcsr: BOOST_LIB_IO=$(BOOST_ROOT)/lib/libboost_iostreams.a
unil-dcsr: BOOST_LIB_PO=$(BOOST_ROOT)/lib/libboost_program_options.a
unil-dcsr: $(BFILE)

static_exe: CXXFLAG=-O2
static_exe: LDFLAG=-O2
static_exe: HTSSRC=../..
static_exe: HTSLIB_INC=$(HTSSRC)/htslib_minimal
static_exe: HTSLIB_LIB=$(HTSSRC)/htslib_minimal/libhts.a
static_exe: BOOST_INC=/usr/include
static_exe: BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
static_exe: BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
static_exe: $(EXEFILE)

#COMPILATION RULES
all: desktop

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

$(EXEFILE): $(OFILE)
	$(CXX) $(LDFLAG) -static -static-libgcc -static-libstdc++ -pthread -o $(EXEFILE) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -Wl,-Bstatic $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE)
