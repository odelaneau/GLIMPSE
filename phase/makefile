#COMPILER MODE C++11
CXX=g++ -std=c++11

#COMPILER & LINKER FLAGS
CXXFLAG=-O3 
LDFLAG=-O3

#DYNAMIC LIBRARIES
DYN_LIBS=-lhts -lboost_iostreams -lboost_program_options -lz -lbz2 -lm -lpthread -llzma

#LCC_phase
BFILE=bin/GLIMPSE_phase
HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#CONDITIONAL PATH DEFINITON

desktop: $(BFILE)

# EXEFILE=bin/GLIMPSE_phase_static
# static_exe: HTSLIB_INC=$(HOME)/Tools/htslib-1.9
# static_exe: HTSLIB_LIB=$(HOME)/Tools/htslib-1.9/libhts.a
# static_exe: BOOST_INC=/usr/include
# static_exe: BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
# static_exe: BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
# static_exe: $(EXEFILE)

#COMPILATION RULES
all: desktop

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^  -o $@ $(DYN_LIBS)

# $(EXEFILE): $(OFILE)
# 	$(CXX) $(LDFLAG) -static -static-libgcc -static-libstdc++ -pthread -o $(EXEFILE) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -Wl,-Bstatic $(DYN_LIBS)
#       $(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc

clean: 
	rm -f obj/*.o $(BFILE)

