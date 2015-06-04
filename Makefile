CC=g++
CFLAGS= -march=native -mtune=native -O3 -Wall -Winline -Wshadow -std=c++11 -fpermissive -fopenmp
#CFLAGS= -fpermissive
LDFLAGS=
SOURCES=mgsolve.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mgsolve
COMMON=

all: clean mgsolve

mgsolve:	
	$(CC) $(CFLAGS) $(SOURCES) -o mgsolve
		
deldata:
	rm -f ./data/*.*
	
clean:
	rm -f *.o mgsolve
	rm -f init.dat
	rm -f solution.dat
		
.PHONY : all clean
