projects = chunk concordance ligate phase split_reference

.PHONY: all $(projects) system clean

all: $(projects)

$(projects):
	$(MAKE) -C $@ $(COMPILATION_ENV)

system:
	for dir in $(projects); do \
	$(MAKE) -C $$dir system; \
	done

clean:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done

