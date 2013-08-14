
#CXX = g++
CXX = g++-mp-4.8
CXXFLAG = -std=c++0x -g
#CXXFLAG = -O3 -std=c++0x -g

#CXX     = /opt/local/bin/clang++-mp-3.3 
#CPPFLAG = -std=c++0x -O3 -g -I${HOME}/include  
#CXXFLAG = -std=c++11 -stdlib=libc++ -I/opt/local/include

LIB_FEMTO  = ~/Drops/femto2
CXXINCLUDE = -I $(LIB_FEMTO)/include
CXXLIB     = -L /opt/local/lib
C_TARGET   = ~/orz/src/sci/ctnew2
F_TARGET   = ${C_TARGET}/fsrc/mrci_tensors

.PHONY: default clean

default: eclipse

SRC = $(notdir $(wildcard *.cc))
OBJ = $(SRC:%.cc=%.o)

eclipse: $(OBJ)
	$(CXX) -o $@ $(OBJ) $(CXXLIB) -L $(LIB_FEMTO) -lfemto -lboost_serialization-mt -lboost_wserialization-mt

%.o: %.cc
	$(CXX) $(CXXFLAG) $(CXXINCLUDE) -c $<

clean: 
	-rm -f *.o
	-rm -f c_*
	-rm -f f_*

install_c:
	-mv c_* ${C_TARGET}

install_f:
	-mv f_* ${C_TARGET}

install:
	-mv c_* ${C_TARGET}
	-mv f_* ${F_TARGET}
