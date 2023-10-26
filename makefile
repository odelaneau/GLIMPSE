projects = chunk concordance ligate phase split_reference extract_num_sites_from_reference_chunk

.PHONY: all $(projects)

all: $(projects)

$(projects):
	$(MAKE) -C $@ $(COMPILATION_ENV)

clean:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done

