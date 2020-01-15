#COMPILER MODE C++11
CXX=g++ -std=c++14

dummy_build_folder_bin := $(shell mkdir -p bin)
dummy_build_folder_obj := $(shell mkdir -p obj)
dummy_build_folder_obj_dbg := $(shell mkdir -p obj/dbg)

#COMPILER & LINKER FLAGS
CXXFLAG=-O3 
LDFLAG=-O3

DBGCXXFLAG=-O0 -g 
DBGLDFLAG=-O0

#DYNAMIC LIBRARIES
DYN_LIBS=-lz -lbz2 -lm -lpthread -llzma

#LCC_phase
BFILE=bin/impute
DBGFILE=bin/impute_dbg
HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)
DBG_OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/dbg/$$(basename $$file .cpp).o; done)

#CONDITIONAL PATH DEFINITON
desktop: HTSLIB_INC=../../htslib
desktop: HTSLIB_LIB=../../htslib/libhts.a
desktop: BOOST_INC=/usr/include
desktop: BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
desktop: BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
desktop: $(BFILE)

cluster: HTSLIB_INC=/software/UHTS/Analysis/samtools/1.4/include
cluster: HTSLIB_LIB=/software/UHTS/Analysis/samtools/1.4/lib64/libhts.a
cluster: BOOST_INC=/software/include
cluster: BOOST_LIB_IO=/software/lib64/libboost_iostreams.a
cluster: BOOST_LIB_PO=/software/lib64/libboost_program_options.a
cluster: $(BFILE)

debug: HTSLIB_INC=../../htslib
debug: HTSLIB_LIB=../../htslib/libhts.a
debug: BOOST_INC=/usr/include
debug: BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
debug: BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
debug: $(DBG_OFILE)
	$(CXX) $(DBGLDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $(DBGFILE) $(DYN_LIBS)

#COMPILATION RULES
all: desktop

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)
	
obj/dbg/%.o: %.cpp $(HFILE)
	$(CXX) $(DBGCXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o obj/dbg/*.o bin/*
