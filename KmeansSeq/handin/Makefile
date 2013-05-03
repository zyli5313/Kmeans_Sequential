TARGET = KmeansSequential
LIBS = -lm
CC = gcc
CFLAGS = -g -Wall

.PHONY: default all clean

default: $(TARGET)

all: default

OBJECTS = ../src/*.c
HEADERS = ../src/*.h

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@


demo:
	../Debug/KmeansSequential -i ../data/input/cluster.csv -o ../data/res/2d.txt -k 2 -t 0.01 -l 20 -d 2
ghc:
	./KmeansSequential -i ../data/input/cluster.csv -o ../data/res/2d.txt -k 2 -t 0.01 -l 20 -d 2

clean:
	-rm -f *.o
	-rm -f $(TARGET)
	
