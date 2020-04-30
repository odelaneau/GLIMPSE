#COMPILER MODE C++11
CXX=g++ -std=c++11

#COMPILER & LINKER FLAGS
CXXFLAG=-O3 
LDFLAG=-O3

#DYNAMIC LIBRARIES
DYN_LIBS=-lhts -lboost_iostreams -lboost_program_options -lz -lbz2 -lm -lpthread -llzma

#LCC_phase
BFILE=bin/GLIMPSE_ligate
HFILE=$(shell find src -name *.h)
CFILE=$(shell find src -name *.cpp)
OFILE=$(shell for file in `find src -name *.cpp`; do echo obj/$$(basename $$file .cpp).o; done)
VPATH=$(shell for file in `find src -name *.cpp`; do echo $$(dirname $$file); done)

#CONDITIONAL PATH DEFINITON
desktop: $(BFILE)

#COMPILATION RULES
all: desktop

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^  -o $@ $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc

clean: 
	rm -f obj/*.o $(BFILE)

