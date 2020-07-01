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
desktop: OPTI_FLAG=-march=native
desktop: HTSLIB_INC=$(HOME)/Tools/htslib-1.9
desktop: HTSLIB_LIB=$(HOME)/Tools/htslib-1.9/libhts.a
desktop: BOOST_INC=/usr/include
desktop: BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
desktop: BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
desktop: $(BFILE)

cluster: OPTI_FLAG=-march=native
cluster: HTSLIB_INC=/software/UHTS/Analysis/samtools/1.4/include
cluster: HTSLIB_LIB=/software/UHTS/Analysis/samtools/1.4/lib64/libhts.a
cluster: BOOST_INC=/software/include
cluster: BOOST_LIB_IO=/software/lib64/libboost_iostreams.a
cluster: BOOST_LIB_PO=/software/lib64/libboost_program_options.a
cluster: $(BFILE)

static_exe: OPTI_FLAG=-mavx2 -mfma
static_exe: HTSLIB_INC=../../htslib
static_exe: HTSLIB_LIB=../../htslib/libhts.a
static_exe: BOOST_INC=/usr/include
static_exe: BOOST_LIB_IO=/usr/local/lib/libboost_iostreams.a
static_exe: BOOST_LIB_PO=/usr/local/lib/libboost_program_options.a
static_exe: $(EXEFILE)

docker: OPTI_FLAG=-mavx2 -mfma
docker: HTSLIB_INC=/home/software/htslib-1.10.2
docker: HTSLIB_LIB=/home/software/htslib-1.10.2/libhts.a
docker: BOOST_INC=/usr/include
docker: BOOST_LIB_IO=/usr/lib/x86_64-linux-gnu/libboost_iostreams.a
docker: BOOST_LIB_PO=/usr/lib/x86_64-linux-gnu/libboost_program_options.a
docker: $(BFILE)

#COMPILATION RULES
all: desktop

$(BFILE): $(OFILE)
	$(CXX) $(LDFLAG) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -o $@ $(DYN_LIBS)

$(EXEFILE): $(OFILE)
	$(CXX) $(LDFLAG) -static -static-libgcc -static-libstdc++ -pthread -o $(EXEFILE) $^ $(HTSLIB_LIB) $(BOOST_LIB_IO) $(BOOST_LIB_PO) -Wl,-Bstatic $(DYN_LIBS)

obj/%.o: %.cpp $(HFILE)
	$(CXX) $(CXXFLAG) $(OPTI_FLAG) -c $< -o $@ -Isrc -I$(HTSLIB_INC) -I$(BOOST_INC)

clean: 
	rm -f obj/*.o $(BFILE)
