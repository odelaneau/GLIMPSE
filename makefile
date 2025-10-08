projects = chunk concordance ligate phase split_reference

BOOST_INC := $(shell brew --prefix boost)/include
BOOST_LIB := $(shell brew --prefix boost)/lib

CXXFLAGS += -I$(BOOST_INC)
LDFLAGS  += -L$(BOOST_LIB)
LDLIBS   += -lboost_program_options

.PHONY: all $(projects) clean

all: $(projects)

$(projects):
	$(MAKE) -C $@ \
		CXXFLAGS="$(CXXFLAGS)" \
		LDFLAGS="$(LDFLAGS)" \
		LDLIBS="$(LDLIBS)"

clean:
	for dir in $(projects); do \
		$(MAKE) -C $$dir clean; \
	done
