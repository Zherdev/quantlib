TARGET = quantlib.a

SRCDIR = src
INCDIR = include
OBJDIR = obj
BINDIR = bin

CC     = mpicc
CFLAGS = -std=gnu99 -Wall -I$(INCDIR)

AR     = ar
ARARGS = rcs

SOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(INCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
rm       = rm -rf

$(BINDIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	@$(AR) $(ARARGS) $(BINDIR)/$(TARGET) $(OBJECTS)
	@echo "Archived "$(BINDIR)/$(TARGET)" successfully!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONY: clean
clean:
	@$(rm) $(BINDIR)
	@$(rm) $(OBJDIR)
	@echo "Cleanup complete!"
