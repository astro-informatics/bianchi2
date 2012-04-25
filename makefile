# Makefile for Bianchi2 library
# Jason McEwen    


# ======== OPTIONS ========

USEPGPLOT = no
#USEPGPLOT = yes


# ======== COMPILER ========

#FC      = nagfor
FC      = ifort

#FC      = gfortran
#FC      = f95
#FC      = g95

ifneq ($(USEPGPLOT),yes)
  OPTPGPLOT     = -DNO_PGPLOT
endif

OPT = $(OPTPGPLOT) $(OPTF95) -DWMAP5 \
      -O3 -DS2_VERSION=\"BIANCHI2_1.0b2\" -DS2_BUILD=\"`svnversion -n .`\" 

OPT = $(OPTPGPLOT) -DMILLIK 

ifeq ($(FC),nagfor)
  OPT += -Wc,-fno-common 
endif
ifeq ($(FC),ifort)
  OPT += -fno-common -openmp
endif
ifeq ($(FC),gfortran)
  OPT += -m64
endif

# NAG ifort Mac OS X libraries typically have 64 bit integers 
# (disable this otherwise).
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  ifeq ($(FC), ifort)
    OPT += -DNAGI8
  endif 
endif


# ======== PPFLAGS ========

ifeq ($(FC),nagfor)
  PPFLAGS = -fpp $(OPT)
else ifeq ($(FC),ifort)
  PPFLAGS = -fpp $(OPT)
else ifeq ($(FC),g95)
  PPFLAGS = -cpp $(OPT)
else ifeq ($(FC),gfortran)
  PPFLAGS = -cpp $(OPT)
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
BIANCHI2SRCPAR= $(BIANCHI2DIR)/src
BIANCHI2SRC   = $(BIANCHI2DIR)/src/mod
BIANCHI2PROG  = $(BIANCHI2DIR)/src/prog
BIANCHI2INC   = $(BIANCHI2DIR)/include
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
    NAGLIBFULL   = /opt/NAG/flmi623dcl/lib/libnag_nag.a \
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


# ======== FFFLAGS ========

FFLAGS  = -I$(HPIXINC) -I$(S2INC) -I$(BIANCHI2INC) -I.


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


# ======== OBJECT FILES TO MAKE ========

BIANCHI2OBJ = $(BIANCHI2INC)/bianchi2_sky_mod.o        \
	      $(BIANCHI2INC)/bianchi2_plm1table_mod.o  \
              $(BIANCHI2INC)/bianchi2_globaldata_mod.o \
              $(BIANCHI2INC)/bianchi2_error_mod.o      


# ======== MAKE RULES ========

default: all

all:     lib prog

lib:	 $(BIANCHI2LIB)/lib$(BIANCHI2LIBNM).a 

prog:    $(BIANCHI2BIN)/bianchi2_sim          \
         $(BIANCHI2BIN)/bianchi2_sim_looptest \
         $(BIANCHI2BIN)/bianchi2_about


$(BIANCHI2INC)/%.o: $(BIANCHI2SRC)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 
	mv *.mod $(BIANCHI2INC)

$(BIANCHI2INC)/%.o: $(BIANCHI2PROG)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 


# Library

$(BIANCHI2LIB)/lib$(BIANCHI2LIBNM).a: $(BIANCHI2OBJ)
	ar -r $(BIANCHI2LIB)/lib$(BIANCHI2LIBNM).a $(BIANCHI2OBJ)


# Documentation

.PHONY: doc
doc:
	doxygen $(BIANCHI2SRCPAR)/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -f $(BIANCHI2DOC)/html/*


# Cleaning up

clean:	tidy
	rm -f $(BIANCHI2INC)/*.mod
	rm -f $(BIANCHI2INC)/*.o
	rm -f $(BIANCHI2LIB)/lib$(BIANCHI2LIBNM).a
	rm -f $(BIANCHI2BIN)/*

tidy:	
	rm -f *.mod
	rm -f $(BIANCHI2SRC)/*~
	rm -f $(BIANCHI2PROG)/*~


# Module dependencies

$(BIANCHI2INC)/bianchi2_error_mod.o:     $(BIANCHI2SRC)/bianchi2_error_mod.f90
$(BIANCHI2INC)/bianchi2_globaldata_mod.o: $(BIANCHI2SRC)/bianchi2_globaldata_mod.f90 \
                                         $(BIANCHI2INC)/bianchi2_error_mod.o
$(BIANCHI2INC)/bianchi2_plm1table_mod.o: $(BIANCHI2SRC)/bianchi2_plm1table_mod.f90 \
                                         $(BIANCHI2INC)/bianchi2_error_mod.o
$(BIANCHI2INC)/bianchi2_sky_mod.o:       $(BIANCHI2SRC)/bianchi2_sky_mod.f90 \
                                         $(BIANCHI2INC)/bianchi2_plm1table_mod.o \
                                         $(BIANCHI2INC)/bianchi2_globaldata_mod.o \
                                         $(BIANCHI2INC)/bianchi2_error_mod.o 


# Program dependencies and compilation

$(BIANCHI2INC)/bianchi2_sim.o: $(BIANCHI2PROG)/bianchi2_sim.f90 lib
$(BIANCHI2BIN)/bianchi2_sim:   $(BIANCHI2INC)/bianchi2_sim.o
	$(FC) -o $(BIANCHI2BIN)/bianchi2_sim $(BIANCHI2INC)/bianchi2_sim.o \
	$(LDFLAGS) $(PPFLAGS)

$(BIANCHI2INC)/bianchi2_sim_looptest.o: $(BIANCHI2PROG)/bianchi2_sim_looptest.f90 lib
$(BIANCHI2BIN)/bianchi2_sim_looptest:   $(BIANCHI2INC)/bianchi2_sim_looptest.o
	$(FC) -o $(BIANCHI2BIN)/bianchi2_sim_looptest $(BIANCHI2INC)/bianchi2_sim_looptest.o \
	$(LDFLAGS) $(PPFLAGS)

$(BIANCHI2ABOUT)/bianchi2_about.o: $(BIANCHI2PROG)/bianchi2_about.f90 lib
	$(FC) -o $(BIANCHI2BIN)/bianchi2_about $(BIANCHI2INC)/bianchi2_about.o \
	$(LDFLAGS) $(PPFLAGS)