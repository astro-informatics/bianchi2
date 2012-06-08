

# ======== OPTIONS ========

USEPGPLOT = no
#USEPGPLOT = yes


# ======== COMPILER ========

#FC      = nagfor
FC      = ifort


ifneq ($(USEPGPLOT),yes)
  OPTPGPLOT     = -DNO_PGPLOT
endif

OPT = $(OPTPGPLOT) -DMILLIK \
      -O3 -DBIANCHI2_VERSION=\"2.0\" -DBIANCHI2_BUILD=\"`svnversion -n .`\" 

ifeq ($(FC),nagfor)
  OPT += -Wc,-fno-common 
endif
ifeq ($(FC),ifort)
  OPT += -fno-common -openmp
endif

# NAG ifort Mac OS X libraries typically have 64 bit integers 
# (disable this otherwise).
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  ifeq ($(FC), ifort)
    OPT += -DNAGI8
  endif 
endif



# ======== LINKS ========

PROGDIR      = ..

HPIXDIR      = $(PROGDIR)/Healpix
ifeq ($(FC),ifort)
  HPIXLIB      = $(HPIXDIR)/lib_ifort
  HPIXINC      = $(HPIXDIR)/include_ifort
else
  ifeq ($(FC),nagfor)
    HPIXLIB      = $(HPIXDIR)/lib_nag
    HPIXINC      = $(HPIXDIR)/include_nag
  else
    HPIXLIB      = $(HPIXDIR)/lib
    HPIXINC      = $(HPIXDIR)/include
  endif
endif
HPIXLIBNM    = healpix

ifeq ($(FC),ifort)
  S2DIR        = $(PROGDIR)/s2_ifort
else
  ifeq ($(FC),nagfor)
    S2DIR        = $(PROGDIR)/s2_nag
  else
    S2DIR        = $(PROGDIR)/s2
  endif
endif
S2LIB        = $(S2DIR)/lib
S2LIBNM      = s2
S2INC        = $(S2DIR)/include
S2SRC        = $(S2DIR)/src/mod
S2PROG       = $(S2DIR)/src/prog
S2BIN        = $(S2DIR)/bin
S2DOC        = $(S2DIR)/doc

BIANCHI2DIR   = $(PROGDIR)/bianchi2
BIANCHI2SRC   = $(BIANCHI2DIR)/src
BIANCHI2MOD   = $(BIANCHI2DIR)/src/mod
BIANCHI2PROG  = $(BIANCHI2DIR)/src/prog
BIANCHI2INC   = $(BIANCHI2DIR)/include
BIANCHI2OBJ   = $(BIANCHI2INC)
BIANCHI2BIN   = $(BIANCHI2DIR)/bin
BIANCHI2LIB   = $(BIANCHI2DIR)/lib
BIANCHI2DOC   = $(BIANCHI2DIR)/doc
BIANCHI2LIBNM = bianchi2

ifeq ($(FC),ifort)
  ifeq ($(UNAME), Darwin)	
    # Mac Air
    CFITSIOLIB   = $(PROGDIR)/cfitsio_ifort/lib    
  else
    # Hypatia
    CFITSIOLIB   = /share/apps/cfitsio/icc/lib
  endif
endif
ifeq ($(FC),nagfor)
  # Mac Air
    CFITSIOLIB   = $(PROGDIR)/cfitsio_nag/lib
endif
CFITSIOLIBNM = cfitsio

ifeq ($(FC),ifort)
  ifeq ($(UNAME), Darwin)	
    # Mac Air

    # Dynamic link.
#    NAGLIBFULL   = -i8 -I/opt/NAG/flmi623dcl/nag_interface_blocks \
#                   /opt/NAG/flmi623dcl/lib/libnag_mkl.a \
#                   -L/opt/NAG/flmi623dcl/mkl_intel64 \
#                   /opt/NAG/flmi623dcl/mkl_intel64/libmkl_intel_ilp64.a \
#                   -lmkl_intel_thread -lmkl_core -liomp5 -lpthread \
#                   -framework IOKit -framework CoreFoundation

    # Static link.
    NAGLIBFULL   = -I/opt/NAG/flmi623dcl/nag_interface_blocks \
                   /opt/NAG/flmi623dcl/lib/libnag_nag.a \
                   -framework IOKit -framework CoreFoundation

  else
    # Hypatia
    NAGLIBFULL   = /share/apps/NAG/fll6i23dcl/lib/libnag_nag.a
  endif
endif
ifeq ($(FC),nagfor)
  # Mac Air
  NAGLIBFULL   = /opt/NAG/flmi622d9l/lib/libnag_nag.a

endif

PGPLOTLIB    = $(PROGDIR)/pgplot
PGPLOTLIBNM  = pgplot
X11LIB       = /usr/X11R6/lib
X11LIBNM     = X11


# ======== SOURCE LOCATIONS ========

vpath %.f90 $(BIANCHI2MOD)
vpath %.f90 $(BIANCHI2PROG)


# ======== FFFLAGS ========

FFLAGS  = -I$(HPIXINC) -I$(S2INC) -I$(BIANCHI2INC) -I.

ifeq ($(FC),nagfor)
  FFLAGS += -mdir $(BIANCHI2INC)
else ifeq ($(FC),ifort)
  FFLAGS += -module $(BIANCHI2INC)
endif


# ======== LDFLAGS ========

ifeq ($(USEPGPLOT),yes)
  LDFLAGSPGPLOT = -L$(PGPLOTLIB) -L$(X11LIB) \
                  -l$(PGPLOTLIBNM) -l$(X11LIBNM)
endif

LDFLAGS =  -L$(BIANCHI2LIB) -l$(BIANCHI2LIBNM) \
           -L$(S2LIB) -l$(S2LIBNM) \
           -L$(HPIXLIB) -l$(HPIXLIBNM) \
           $(NAGLIBFULL) \
           -L$(CFITSIOLIB) -l$(CFITSIOLIBNM) \
           $(LDFLAGSPGPLOT)
#           -L$(NAGLIB) -l$(NAGLIBNM) \


# ======== PPFLAGS ========

ifeq ($(FC),nagfor)
  PPFLAGS = -fpp $(OPT)
else ifeq ($(FC),ifort)
  PPFLAGS = -fpp $(OPT)
endif


# ======== OBJECT FILES TO MAKE ========

BIANCHI2OBJS = $(BIANCHI2OBJ)/bianchi2_sky_mod.o        \
	       $(BIANCHI2OBJ)/bianchi2_lut_mod.o        \
               $(BIANCHI2OBJ)/bianchi2_globaldata_mod.o \
               $(BIANCHI2OBJ)/bianchi2_error_mod.o      

BIANCHI2BINS = $(BIANCHI2BIN)/bianchi2_about            \
               $(BIANCHI2BIN)/bianchi2_sim              \
               $(BIANCHI2BIN)/bianchi2_lut_gen


# ======== MAKE RULES ========

$(BIANCHI2OBJ)/%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(PPFLAGS) 

$(BIANCHI2BIN)/%: %.f90 $(BIANCHI2LIB)/lib$(BIANCHI2LIBNM).a
	$(FC) $< -o $@ $(FFLAGS) $(LDFLAGS) $(PPFLAGS)  


# Default

.PHONY: default
default: all

.PHONY: all
all:     lib prog


# Library

.PHONY: lib
lib:	 $(BIANCHI2LIB)/lib$(BIANCHI2LIBNM).a
$(BIANCHI2LIB)/lib$(BIANCHI2LIBNM).a: $(BIANCHI2OBJS)
	ar -r $(BIANCHI2LIB)/lib$(BIANCHI2LIBNM).a $(BIANCHI2OBJS)


# Programs

.PHONY: prog
prog: $(BIANCHI2BINS)


# Tests

.PHONY: runtest
runtest: prog
	./bin/bianchi2_sim param.par


# Documentation

.PHONY: doc
doc:
	doxygen $(BIANCHI2SRC)/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -f $(BIANCHI2DOC)/html/*


# Cleaning up

.PHONY: clean
clean:	tidy
	rm -f $(BIANCHI2INC)/*.mod
	rm -f $(BIANCHI2INC)/*.o
	rm -f $(BIANCHI2LIB)/lib$(BIANCHI2LIBNM).a
	rm -f $(BIANCHI2BIN)/*

.PHONY: tidy
tidy:	
	rm -f *.mod
	rm -f $(BIANCHI2MOD)/*~
	rm -f $(BIANCHI2PROG)/*~


# Module dependencies

$(BIANCHI2OBJ)/bianchi2_error_mod.o:      $(BIANCHI2MOD)/bianchi2_error_mod.f90
$(BIANCHI2OBJ)/bianchi2_globaldata_mod.o: $(BIANCHI2MOD)/bianchi2_globaldata_mod.f90 \
                                          $(BIANCHI2OBJ)/bianchi2_error_mod.o
$(BIANCHI2OBJ)/bianchi2_lut_mod.o:        $(BIANCHI2MOD)/bianchi2_lut_mod.f90        \
                                          $(BIANCHI2OBJ)/bianchi2_error_mod.o
$(BIANCHI2OBJ)/bianchi2_sky_mod.o:        $(BIANCHI2MOD)/bianchi2_sky_mod.f90        \
                                          $(BIANCHI2OBJ)/bianchi2_lut_mod.o          \
                                          $(BIANCHI2OBJ)/bianchi2_globaldata_mod.o   \
                                          $(BIANCHI2OBJ)/bianchi2_error_mod.o 

