projects = chunk concordance inspect ligate phase split_reference

.PHONY: all $(projects) system clean

all: $(projects)

$(projects):
	$(MAKE) -C $@ $(COMPILATION_ENV)

system: $(addprefix system-,$(projects))
system-%:
	$(MAKE) -C $* system

clean: $(addprefix clean-,$(projects))
clean-%:
	$(MAKE) -C $* clean

