# ------------------------------------------------------------
# Makefile for the unified 2D UGKS + Hermite-spectral solver
# Supported compilers: gfortran (default) and Intel ifort.
# Select Intel Fortran with: make FC=ifort
#
# Fortran modules require ordered compilation, so this Makefile compiles
# the sources in a fixed dependency order via a single serial rule.
# ------------------------------------------------------------
FC = gfortran

SRCDIR = src
EXEC = ugks2d

FC_NAME := $(notdir $(firstword $(FC)))
ifneq ($(filter ifort ifort.exe,$(FC_NAME)),)
FFLAGS ?= -O2
OMPFLAG ?= -qopenmp
MODFLAGS = -module $(SRCDIR) -I$(SRCDIR)
else
FFLAGS ?= -O2 -ffree-line-length-none -fno-range-check
OMPFLAG ?= -fopenmp
MODFLAGS = -J$(SRCDIR) -I$(SRCDIR)
endif

# Compile order = module dependency order (each USEs the ones above it).
SRCS = config data2d grid gastheory2d gravity hermite_symmetric2d velocity boundary flux \
       initialization solver output2d main
OBJS = $(addprefix $(SRCDIR)/,$(addsuffix .o,$(SRCS)))

.PHONY: all clean run run5 run6

all: $(EXEC)

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) $(OMPFLAG) $(OBJS) -o $(EXEC)

# Serial compile of all sources in dependency order.
$(OBJS): $(SRCDIR)/.compile.stamp ;

$(SRCDIR)/.compile.stamp: $(addprefix $(SRCDIR)/,$(addsuffix .f90,$(SRCS)))
	@rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod
	@for f in $(SRCS); do \
	  echo "  FC  $(SRCDIR)/$$f.f90"; \
	  $(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/$$f.f90 -o $(SRCDIR)/$$f.o || exit 1; \
	done
	@touch $@

run: run6
run5: $(EXEC)
	./$(EXEC) cases/poiseuille_fast.namelist
run6: $(EXEC)
	./$(EXEC) cases/rt_60_fast.namelist

clean:
	-rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/.compile.stamp $(EXEC)
# ------------------------------------------------------------
