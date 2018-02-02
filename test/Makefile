# (C) Copyright Kirill Lykov 2013.
#
# Distributed under the FreeBSD Software License (See accompanying file license.txt) 
#
# Make file for oni2avi. Check include and lib pathes before runing.

SHELL = /bin/sh
.SUFFIXES: .cpp .u
.IGNORE:

ROOT =  runtest
EXE =   $(ROOT)

UNAME := $(shell uname)
CC = g++

INCLUDES = -I/home/klykov/bin/gtest/include -I/home/klykov/bin/hdf5/include -I../include
LIBS_PATH = -L/home/klykov/bin/gtest/lib -L/home/klykov/bin/hdf5/lib

LIBS = -lgtest -lpthread -lhdf5 -lhdf5_cpp

CCFLAGS = -std=c++11 -O0 -g -DUSE_HDF5
LINK = g++
SRC = test_main.cpp
OBJ = test_main.o

$(EXE): $(OBJ)
	$(LINK) -o $(EXE) $(OBJ) $(LIBS_PATH) $(LIBS)

%.o:%.cpp
	$(CC) $(CCFLAGS) $(INCLUDES) -c $<

clean:
	rm -fr $(ROOT) *.o *.d *.dSYM
