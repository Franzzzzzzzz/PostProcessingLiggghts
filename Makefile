CC=g++

# Basic compilation
#CFLAGS=-c -O3 -I./Headers -I/opt/local/include
#LDFLAGS=-L./newmat10 
#LDFLAGS=-L/opt/local/lib
#LIBS=-lm -lpthread 
#-lnewmatMac

# Compilation incluant le support matlab
#CFLAGS=-c -O3 -m32 -DMATLAB -I./Headers -I/opt/local/include -I/home/franz/Logiciels/matlab2009/extern/include
CFLAGS=-c -g -DMATLAB -DBOOST -DUSETIFF -I./Headers -I/opt/local/include -I/Applications/MATLAB_R2015b.app/extern/include
# -DUSETIFF and -ltiff to use tiff
#LDFLAGS=-L/opt/local/lib -L/home/franz/Logiciels/matlab2009/bin/glnx86
LDFLAGS=-L/opt/local/lib -L/Applications/MATLAB_R2015b.app/bin/maci64
LIBS=-lstdc++ -lm -lpthread -lmat -lmx -lz -ltiff

SOURCES=Ids.cpp Gunzip.cpp Writing.cpp Calculs.cpp Coarse.cpp Compress.cpp Dump.cpp FDump.cpp LDump.cpp LcfDump.cpp Filter.cpp Main.cpp Statistics.cpp Step.cpp FStep.cpp LStep.cpp CFStep.cpp Surface.cpp Multisphere.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=PostProcessing

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean : 
	rm *.o
