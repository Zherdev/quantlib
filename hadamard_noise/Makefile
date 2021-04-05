TARGET   = hadamard

SRCDIR    = src
INCDIR    = include
OBJDIR    = obj
BINDIR    = bin
DATADIR   = data
REPORTDIR = report

CC       = mpixlc_r
CFLAGS   = -qsmp=omp -lm -I$(INCDIR) -O2

LINKER   = mpixlc_r
LFLAGS   = -qsmp=omp -lm -I$(INCDIR)

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

QUBITS := 24
.PHONY: generate
generate: $(BINDIR)/$(TARGET)
	@mkdir -p $(DATADIR)
	for Q in $(QUBITS); do \
		mpisubmit.bg -n 16 "bin/hadamard" -- generate $$Q 4 $(DATADIR)/input-$$Q.bin ; \
	done
	@echo "Generation started!"

.PHONY: time
time: $(BINDIR)/$(TARGET)
	bash bluegene_time.sh
	@echo "Time run is done!"

.PHONY: loss
loss: $(BINDIR)/$(TARGET)
	bash bluegene_loss.sh
	@echo "Loss run is done!"
