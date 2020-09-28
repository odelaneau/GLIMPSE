#COMPILER MODE C++11
CXX=g++ -std=c++11

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

cluster: HTSSRC=/software/UHTS/Analysis/samtools/1.4
cluster: HTSLIB_INC=$(HTSSRC)/include
cluster: HTSLIB_LIB=$(HTSSRC)/lib64/libhts.a
cluster: BOOST_INC=/software/include
cluster: BOOST_LIB_IO=/software/lib64/libboost_iostreams.a
cluster: BOOST_LIB_PO=/software/lib64/libboost_program_options.a
cluster: $(BFILE)

static_exe: HTSSRC=../..
static_exe: HTSLIB_INC=$(HTSSRC)/htslib
static_exe: HTSLIB_LIB=$(HTSSRC)/htslib/libhts.a
static_exe: BOOST_INC=/usr/include
static_exe: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
static_exe: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
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
