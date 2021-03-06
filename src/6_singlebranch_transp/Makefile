# ====================================================================
#   "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2014-2015
#
#                    Copyright D. Notaro 2015
# ====================================================================
#   FILE        : Makefile
#   DESCRIPTION : makefile for test simulations
#   AUTHOR      : Domenico Notaro <domenico.not@gmail.com>
#   DATE        : April 2015
# ====================================================================

CPPFLAGS=-I../../include -I$(mkGetfemInc)
CXXFLAGS=-std=c++11 -D=M3D1D_VERBOSE_
CXXFLAGS += -I ${SAMG}/
CXXFLAGS+= -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE -DPYRAMID_TRIANGULAR_FACETS
ifeq ($(DEBUG),yes)
  OPTFLAGS=-g -Wall
else
  OPTFLAGS=-O3 -march=native
  CPPFLAGS+=-DNDEBUG
endif
LDFLAGS=-L../../lib -L$(mkGetfemLib)
LDFLAGS += -L/opt/lib/samg/
LIBRARIES=-lproblem3d1d -lgetfem
LIBRARIES += -lamg -liomp5
LINKLIBS := ../../lib/libproblem3d1d.a


SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
EXEC=M3D1D

OUTDIR=vtk

.PHONY: all clean distclean

all: library $(EXEC) 
	@echo
	@echo Compilation completed!

%.o: %.cpp $(LINKLIBS) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -o $@ -c $<

$(EXEC): $(OBJS) 
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBRARIES)

library:
	$(MAKE) -C ../../include

clean:
	$(RM) $(OBJS) $(EXEC) *~ *.log

distclean: clean
	$(RM) *.txt $(OUTDIR)/*
