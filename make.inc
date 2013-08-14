
SHELL = /bin/sh

# Home directory for FEMTO
FEMTO_HOME=~/Drops/femto2

# Include sirectories
CPPINCLUDE = -I $(FEMTO_HOME)/include

# C++ compiler specification
#CPP      = g++
CPP      = g++-mp-4.8
#CPPFLAG  = -O3 -std=c++0x -g
CPPFLAG = -std=c++0x -g

#CPP     = /opt/local/bin/clang++-mp-3.3 
#CPPFLAG = -std=c++0x -O3 -g -I${HOME}/include  
#CPPFLAG = -std=c++11 -stdlib=libc++ -I/opt/local/include

# The archiver specification
ARCH      = ar
APCHFLAGS = cr
RANLIB    = ranlib

# Name of the resultant FEMTO library
FEMTOLIB = libfemto.a
