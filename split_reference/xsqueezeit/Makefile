# TODO : Set this path accordingly
ifeq ($(HTSLIB_LIB),)
HTSLIB_LIB := /home/srubinac/git/htslib
endif
ifeq ($(HTSLIB_INC),)
HTSLIB_INC := /home/srubinac/git/htslib/htslib
endif
ifeq ($(ZSTD_LIB),)
ZSTD_LIB := ../zstd/lib
endif
ifeq ($(ZSTD_INC),)
ZSTD_INC := ../zstd/lib
endif

# C++ Compiler
CXX=g++
INCLUDE_DIRS=-I include -I $(HTSLIB_INC) -I $(ZSTD_INC)
#EXTRA_FLAGS=-fsanitize=address -fsanitize=undefined -fsanitize=pointer-subtract -fsanitize=pointer-compare -fno-omit-frame-pointer -fstack-protector-all -fcf-protection
CXXFLAGS=-O3 -g -Wall -std=c++17 $(INCLUDE_DIRS) $(CXXEXTRAFLAGS) $(EXTRA_FLAGS)
# Linker
LD=g++
LIBS=-lpthread -lhts -lzstd
LDFLAGS=-O3 $(EXTRA_FLAGS) -L $(HTSLIB_LIB) -L $(ZSTD_LIB)

# Project specific :
CPP_SOURCES := $(wildcard *.cpp)
CPP_OBJS := $(CPP_SOURCES:.cpp=.o)
DEPENDENCIES := $(CPP_SOURCES:.cpp=.d)
DEPENDENCIES := $(DEPENDENCIES:.c=.d)

# Rules
all : $(CPP_OBJS)
	ar -rsv libxsqueezeit.a $^

# Do not include the depency rules for "clean"
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPENDENCIES)
endif

# Compile
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to generate the dependency files
%.d : %.cpp
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@
%.d : %.c
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@

# Remove artifacts
clean :
	rm -f $(CPP_OBJS) $(DEPENDENCIES)
	rm -f libxsqueezeit.a

# Rules that don't generate artifacts
.PHONY : all clean
