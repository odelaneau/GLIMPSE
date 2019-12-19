#module load Bioinformatics/Software/vital-it
#module add Development/gcc/9.2.1

#COMPILER MODE C++11
CXX=g++ -std=c++0x

#HTSLIB LIBRARY [SPECIFY YOUR OWN PATHS]
HTSLIB_INC=$(HOME)/Tools/htslib-1.9
HTSLIB_LIB=$(HOME)/Tools/htslib-1.9/libhts.a

#BOOST IOSTREAM & PROGRAM_OPTION LIBRARIES [SPECIFY YOUR OWN PATHS]
BOOST_INC=/usr/include
BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a

#HTSLIB LIBRARY [SPECIFY YOUR OWN PATHS]
#HTSLIB_INC=/software/UHTS/Analysis/samtools/1.4/include
#HTSLIB_LIB=/software/UHTS/Analysis/samtools/1.4/lib64/libhts.a

#BOOST IOSTREAM & PROGRAM_OPTION LIBRARIES [SPECIFY YOUR OWN PATHS]
#BOOST_INC=/software/include
#BOOST_LIB_IO=/software/lib64/libboost_iostreams.a
#BOOST_LIB_PO=/software/lib64/libboost_program_options.a



#COMPILER & LINKER FLAGS
CXXFLAG=-O3
LDFLAG=-O3
#CXXFLAG=-g
#LDFLAG=-g

#DYNAMIC LIBRARIES
DYN_LIBS=-lz -lbz2 -lm -lpthread -llzma

#SHAPEIT SOURCES & BINARY
BFILE=bin/LCC-phase
HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#COMPILATION RULES
all: $(BFILE)

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE)
