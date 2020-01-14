#COMPILER MODE C++11
CXX=g++ -std=c++11

#COMPILER & LINKER FLAGS
CXXFLAG=-O3 
LDFLAG=-O3

#DYNAMIC LIBRARIES
DYN_LIBS=-lz -lbz2 -lm -lpthread -llzma

#LCC_phase
BFILE=bin/LCC_ligate
HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#CONDITIONAL PATH DEFINITON
desktop: HTSLIB_INC=$(HOME)/Tools/htslib-1.9
desktop: HTSLIB_LIB=$(HOME)/Tools/htslib-1.9/libhts.a
desktop: BOOST_INC=/usr/include
desktop: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
desktop: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
desktop: $(BFILE)

axiom: module load Bioinformatics/Software/vital-it
axiom: module load module add Development/gcc/9.2.1
axiom: HTSLIB_INC=/software/UHTS/Analysis/samtools/1.4/include
axiom: HTSLIB_LIB=/software/UHTS/Analysis/samtools/1.4/lib64/libhts.a
axiom: BOOST_INC=/software/include
axiom: BOOST_LIB_IO=/software/lib64/libboost_iostreams.a
axiom: BOOST_LIB_PO=/software/lib64/libboost_program_options.a
axiom: $(BFILE)

jura: HTSLIB_INC=/software/UHTS/Analysis/samtools/1.4/include
jura: HTSLIB_LIB=/software/UHTS/Analysis/samtools/1.4/lib64/libhts.a
jura: BOOST_INC=/software/include
jura: BOOST_LIB_IO=/software/lib64/libboost_iostreams.a
jura: BOOST_LIB_PO=/software/lib64/libboost_program_options.a
jura: $(BFILE)

#COMPILATION RULES
all: desktop

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE)
