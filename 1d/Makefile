# ------------------------------------------------------------
# Makefile for the unified 1D UGKS + Hermite-spectral solver
# Supported compilers: gfortran (default) and Intel ifort.
# Select Intel Fortran with: make FC=ifort
#
# Fortran modules require ordered compilation, so this Makefile compiles
# the sources in a fixed dependency order via a single serial rule.
# ------------------------------------------------------------
FC = gfortran

SRCDIR = src
EXEC_BASE = ugks1d

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
SRCS = config data1d grid gastheory gravity velocity flux initialization solver output main
OBJS = $(addprefix $(SRCDIR)/,$(addsuffix .o,$(SRCS)))

.PHONY: all clean run run2 run3 run4

ifeq ($(OS),Windows_NT)
EXEC = $(EXEC_BASE).exe
CLEAN_BUILD = powershell.exe -NoProfile -Command "Remove-Item -Path '$(SRCDIR)/*.o','$(SRCDIR)/*.mod','$(SRCDIR)/.compile.stamp','$(EXEC)' -Force -ErrorAction SilentlyContinue"
TOUCH_STAMP = powershell.exe -NoProfile -Command "Set-Content -NoNewline -Path '$(SRCDIR)/.compile.stamp' -Value ''"
else
EXEC = $(EXEC_BASE)
CLEAN_BUILD = rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/.compile.stamp $(EXEC)
TOUCH_STAMP = touch $(SRCDIR)/.compile.stamp
endif

all: $(EXEC)

# Link only (objects built by the serial rule below)
$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) $(OMPFLAG) $(OBJS) -o $(EXEC)

# Serial compile of all sources in dependency order.
# Touching any source rebuilds from that point onward (stamp file tracks it).
$(OBJS): $(SRCDIR)/.compile.stamp ;

$(SRCDIR)/.compile.stamp: $(addprefix $(SRCDIR)/,$(addsuffix .f90,$(SRCS)))
	$(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/config.f90 -o $(SRCDIR)/config.o
	$(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/data1d.f90 -o $(SRCDIR)/data1d.o
	$(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/grid.f90 -o $(SRCDIR)/grid.o
	$(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/gastheory.f90 -o $(SRCDIR)/gastheory.o
	$(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/gravity.f90 -o $(SRCDIR)/gravity.o
	$(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/velocity.f90 -o $(SRCDIR)/velocity.o
	$(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/flux.f90 -o $(SRCDIR)/flux.o
	$(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/initialization.f90 -o $(SRCDIR)/initialization.o
	$(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/solver.f90 -o $(SRCDIR)/solver.o
	$(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/output.f90 -o $(SRCDIR)/output.o
	$(FC) $(FFLAGS) $(OMPFLAG) $(MODFLAGS) -c $(SRCDIR)/main.f90 -o $(SRCDIR)/main.o
	$(TOUCH_STAMP)

run: run3
run2: $(EXEC)
	./$(EXEC) cases/hydrostatic_kn1_gh36.namelist
run3: $(EXEC)
	./$(EXEC) cases/sod_kn001_gh34.namelist
run4: $(EXEC)
	./$(EXEC) cases/fourier_force2_kn01_gh22.namelist

clean:
	-$(CLEAN_BUILD)
# ------------------------------------------------------------
