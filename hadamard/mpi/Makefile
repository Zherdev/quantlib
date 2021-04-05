TARGET   = hadamard

SRCDIR    = src
INCDIR    = include
OBJDIR    = obj
BINDIR    = bin
DATADIR   = data
REPORTDIR = report

CC       = mpicc
CFLAGS   = -std=gnu99 -Wall -I$(INCDIR) -lm

LINKER   = mpicc
LFLAGS   = -Wall -I$(INCDIR) -lm

SOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(INCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
rm       = rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	@$(LINKER) $(OBJECTS) $(LFLAGS) -o $@
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONY: clean
clean:
	@$(rm) $(BINDIR)/$(TARGET)
	@$(rm) $(OBJECTS)
	@$(rm) -f $(REPORTDIR)/*
	@$(rm) *.out *.err
	@echo "Cleanup complete!"

Ns := 25 26 27 28
.PHONY: generate
generate: $(BINDIR)/$(TARGET)
	@mkdir -p $(DATADIR)
	for N in $(Ns); do \
		mpisubmit-new.pl -p 16 "bin/hadamard" -- generate $$N $(DATADIR)/input-$$N.bin ; \
	done
	@echo "Generation started!"

.PHONY: test
test: $(BINDIR)/$(TARGET)
	bash polus_test.sh
	@echo "Test is done!"

.PHONY: run
run: $(BINDIR)/$(TARGET)
	bash polus_run.sh
	@echo "Run is done!"
