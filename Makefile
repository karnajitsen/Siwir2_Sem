CC=g++
CFLAGS= -O3 -Wall -Winline -Wshadow -std=c++11 -fpermissive
#CFLAGS= -fpermissive
LDFLAGS=
SOURCES=mgsolve.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mgsolve
COMMON=

all: clean mgsolve

mgsolve:
	$(CC) $(CFLAGS) $(SOURCES) -o mgsolve
	
test: clean
	./mgsolve 8 5
	git add data
	git commit -m "new data file"

onegrid: deldata
	./mgsolve 8 10

allgrid: deldata
	./mgsolve 3 20
	./mgsolve 4 20
	./mgsolve 5 20
	./mgsolve 6 20
	./mgsolve 7 20
	./mgsolve 8 20
	
deldata:
	rm -f ./data/Neumann/*.txt
	rm -f ./data/Dirichlet/*.txt

clean:
	rm -f *.o mgsolve
	rm -f ./data/Neumann/*.txt
	rm -f ./data/Dirichlet/*.txt
	
.PHONY : all clean
