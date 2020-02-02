projects = chunk impute ligate concordance sample

.PHONY: all $(projects)

all: $(projects)

$(projects):
	$(MAKE) -C $@

clean:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done

