ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
ROOTLIBDIR = $(shell root-config --libdir)

GXX      = g++ -Wall -O3 -DNDEBUG


LINKROOT = $(ROOTCFLAGS) -L$(ROOTLIBDIR) $(ROOTGLIBS) $(ROOTLIBS) -lASImage

all: simpix_start simpix

simpix_start: simpix_start.cpp
	$(GXX) -o simpix_start simpix_start.cpp $(LINKROOT)

simpix: simpix.cpp
	$(GXX) -o simpix simpix.cpp $(LINKROOT)

clean:
	rm -f simpix_start simpix out.png newout.png collage.png newcollage.png

