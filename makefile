projects = chunk concordance ligate phase sample snparray

.PHONY: all $(projects)

all: $(projects)

$(projects):
	$(MAKE) -C $@ $(COMPILATION_ENV)

clean:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done

