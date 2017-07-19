CC=g++

# Basic compilation
#CFLAGS=-c -O3 -I./Headers -I/opt/local/include
#LDFLAGS=-L./newmat10 
#LDFLAGS=-L/opt/local/lib
#LIBS=-lm -lpthread 

# Compilation incluant le support matlab
#CFLAGS=-c -O3 -m32 -DMATLAB -I./Headers -I/opt/local/include -I/home/franz/Logiciels/matlab2009/extern/include
# -DUSETIFF and -ltiff to use tiff
#LDFLAGS=-L/opt/local/lib -L/home/franz/Logiciels/matlab2009/bin/glnx86


SOURCES=Ids.cpp Gunzip.cpp Writing.cpp Calculs.cpp Coarse.cpp Compress.cpp Dump.cpp FDump.cpp LDump.cpp LcfDump.cpp Filter.cpp Main.cpp Statistics.cpp Step.cpp FStep.cpp LStep.cpp CFStep.cpp Surface.cpp Multisphere.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=PostProcessing

all: CFLAGS=-c -O3 -I./Headers 
all: LDFLAGS=
all: LIBS=-lstdc++ -lm -lpthread -lz
all: info $(SOURCES) $(EXECUTABLE)


advanced: CFLAGS=-c -O3 -DMATLAB -DBOOST -DUSETIFF -I./Headers -I/opt/local/include -I/Applications/MATLAB_R2015b.app/extern/include
advanced: LDFLAGS=-L/opt/local/lib -L/Applications/MATLAB_R2015b.app/bin/maci64
advanced: LIBS=-lstdc++ -lm -lpthread -lmat -lmx -lz -ltiff
advanced: $(SOURCES) $(EXECUTABLE)

debug: CFLAGS=-c -g -DMATLAB -DBOOST -DUSETIFF -I./Headers -I/opt/local/include -I/Applications/MATLAB_R2015b.app/extern/include
debug: LDFLAGS=-L/opt/local/lib -L/Applications/MATLAB_R2015b.app/bin/maci64
debug: LIBS=-lstdc++ -lm -lpthread -lmat -lmx -lz -ltiff
debug: $(SOURCES) $(EXECUTABLE)

info:
	@echo 'Compilation without Matlab or Boost or Tiff support. Use "make advanced" to activate them.'

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean : 
	rm *.o
