
# project name and executable name
project = calib
# source suffix
source = .cpp 
# object files to make
objects = Engine.o TofCalibrationMaker.o PicoSplitter.o SplineMaker.o TofTuple.o

sources = $(wildcard *.cpp)
headers = $(wildcard ../include/*.h)

# ROOT libs and includes
ROOTCFLAGS    	= $(shell root-config --cflags)
ROOTLDFLAGS    	= $(shell root-config --ldflags)
ROOTLIBS      	= $(shell root-config --libs)
ROOTGLIBS     	= $(shell root-config --glibs)

# Custom Project includes and libs
# -lGui -lHtml -lGuiHtml
includes		= -I../include -I$(JDB_LIB)/include/jdb -I$(JDB_LIB)/include/ -I. -I$(ROOTSYS)/include
libraries 		= $(ROOTLIBS) $(SYSLIBS) -lMathCore -lMathMore
staticLibs 		= $(JDB_LIB)/lib/libJDB.a

cxx 		= g++ -std=c++0x
flags 		= -Wall -O3 $(ROOTCFLAGS) $(includes) -Wno-write-strings
compile 	= $(cxx) $(flags) -c 
ldFlags  	= $(ROOTLDFLAGS)

all: $(project)

$(project): $(objects)
		$(cxx) $(ldFlags) -o  $(project) $(objects) $(staticLibs) $(ROOTGLIBS) $(libraries)
		@cp $(project) ../bin

%.o: %$(source)
		$(compile)  $<

clean:
		@rm -f $(objects) $(project)
		@rm -f ../bin/$(project)
		@echo "clean done"

full:
	clean
	$(project)
