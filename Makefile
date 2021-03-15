CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.cpp output.cpp Born.cpp delta.cpp general.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=start

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -rf *.o *.csv $(EXECUTABLE)