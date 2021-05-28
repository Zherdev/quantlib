TARGET = cnot

SOURCE   = cnot.c
OBJECT   = cnot.o
LIBDIR   = ../../lib
INCDIR   := $(LIBDIR)/include
QUANTLIB := $(wildcard $(LIBDIR)/bin/*.a)

CC     = mpicc
CFLAGS = -std=gnu99 -Wall -I$(INCDIR)
rm     = rm -rf

LINKER   = mpicc
LFLAGS   = -Wall -I$(INCDIR) -lm

$(TARGET): $(OBJECT)
	@$(LINKER) $(OBJECT) $(QUANTLIB) $(LFLAGS) -o $@
	@echo "Linking complete!"

$(OBJECT): $(SOURCE) $(QUANTLIB)
	@$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONY: clean
clean:
	@$(rm) $(OBJECT)
	@$(rm) $(TARGET)
	@echo "Cleanup complete!"
