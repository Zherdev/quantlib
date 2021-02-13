TARGET   = adamar

SRCDIR   = src
INCDIR   = include
OBJDIR   = obj
BINDIR   = bin

CC       = gcc
CFLAGS   = -std=gnu99 -Wall -I$(INCDIR) -fopenmp -lm

LINKER   = gcc
LFLAGS   = -Wall -I$(INCDIR) -fopenmp -lm

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
	@$(rm) -f report/*
	@$(rm) *.out *.err
	@echo "Cleanup complete!"