MAIN_FILE=main
LIB_PATH=../../lib/
LIB_NAME=matrix-lib

all: main

main: $(MAIN_FILE).o
	g++ $(MAIN_FILE).o -o main

$(MAIN_FILE).o: $(MAIN_FILE).cpp
	g++ -c $(MAIN_FILE).cpp

clean:
	rm -rf *o main