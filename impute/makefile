#COMPILER MODE C++11
CXX=g++ -std=c++0x

dummy_build_folder_bin := $(shell mkdir -p bin)
dummy_build_folder_obj := $(shell mkdir -p obj)
dummy_build_folder_obj_dbg := $(shell mkdir -p obj/dbg)

#HTSLIB LIBRARY [SPECIFY YOUR OWN PATHS]
HTSLIB_INC=../htslib
HTSLIB_LIB=../htslib/libhts.a

#BOOST IOSTREAM & PROGRAM_OPTION LIBRARIES [SPECIFY YOUR OWN PATHS]
BOOST_INC=/usr/local/include/boost
BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a

#COMPILER & LINKER FLAGS
CXXFLAG=-O3
LDFLAG=-O3
DBGCXXFLAG=-g
DBGLDFLAG=-g

#DYNAMIC LIBRARIES
DYN_LIBS=-lz -lbz2 -lm -lpthread -llzma -lcurl -lssl -lcrypto #-lz -lbz2 -lm -lpthread -llzma

#SOURCES & BINARY
BFILE=bin/LCC_impute
EXE=bin/LCC_impute_static
DBGFILE=bin/LCC_impute_debug

HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)
DBG_OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/dbg/$$(basename $$file .cpp).o; done)

#COMPILATION RULES
all: opti

opti: $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $(BFILE) $(DYN_LIBS)

#exe requires htslib's library to be compiled in a light mode, disabling libcurl, etc.:
#autoheader && autoconf && ./configure --disable-lzma --disable-bz2 --disable-libcurl && make libhts.a
exe: $(OFILE)
	$(CXX) $(LDFLAG) -static -static-libgcc -static-libstdc++ -pthread -o $(EXE) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -Wl,-Bstatic $(DYN_LIBS)

debug: $(DBG_OFILE)
	$(CXX) $(DBGLDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $(DBGFILE) $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

obj/dbg/%.o: %.cpp $(HFILE)
	$(CXX) $(DBGCXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o obj/dbg/*.o bin/*
