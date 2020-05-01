# projects = chunk concordance impute ligate phase sample snparray
projects = chunk ligate phase sample

.PHONY: all $(projects)

all: $(projects)

$(projects):
	$(MAKE) -C $@

clean:
	for dir in $(projects); do \
	$(MAKE) $@ -C $$dir; \
	done

