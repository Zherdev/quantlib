TARGET = test

SRCDIR    = src
INCDIR    = include
LIBDIR    = ../lib
OBJDIR    = obj
BINDIR    = bin
REPDIR    = report

SOURCES   := $(wildcard $(SRCDIR)/*.c)
INCLUDES  := $(wildcard $(INCDIR)/*.h)
QUANTLIB  := $(wildcard $(LIBDIR)/$(BINDIR)/*.a)
OBJECTS   := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
rm        = rm -rf
MPIEXEC   = mpiexec
MPISUBMIT = mpisubmit.pl

CC       = mpicc
CFLAGS   = -std=gnu99 -Wall -I$(INCDIR) -I$(LIBDIR)/$(INCDIR)

LINKER   = mpicc
LFLAGS   = -Wall -I$(INCDIR) -I$(LIBDIR)/$(INCDIR) -lm

$(BINDIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	@$(LINKER) $(OBJECTS) $(QUANTLIB) $(LFLAGS) -o $@
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c $(QUANTLIB)
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONY: clean
clean:
	@$(rm) $(BINDIR)
	@$(rm) $(OBJDIR)
	@$(rm) $(REPDIR)
	@echo "Cleanup complete!"

.PHONY: run-mpiexec
run-mpiexec: $(BINDIR)/$(TARGET)
	for procs in 1 2 4 ; do                                        \
		echo "";                                                   \
		echo "Run tests on $$procs processes";                     \
		echo "";                                                   \
    	$(MPIEXEC) -n $$procs --oversubscribe $(BINDIR)/$(TARGET); \
	done

.PHONY: run-mpisubmit
run-mpisubmit: $(BINDIR)/$(TARGET)
	@mkdir -p $(REPDIR)
	for procs in 1 2 4 ; do                    \
		echo "";                               \
		echo "Run tests on $$procs processes"; \
		echo "";                               \
    	$(MPISUBMIT) -p $$procs --stdout $(REPDIR)/$$procs.out --stderr $(REPDIR)/$$procs.err  $(BINDIR)/$(TARGET); \
		echo "See quantlib/test/$(REPDIR) for test results."; \
	done
