MAIN_FILE=main
LIB_PATH=../../lib/
LIB_NAME=matrix-lib

all: main

main: $(MAIN_FILE).o $(LIB_NAME).o
	g++ $(MAIN_FILE).o $(LIB_NAME).o -o main

$(MAIN_FILE).o: $(MAIN_FILE).cpp
	g++ -c $(MAIN_FILE).cpp

$(LIB_NAME).o: $(LIB_PATH)$(LIB_NAME).cpp
	g++ -c $(LIB_PATH)$(LIB_NAME).cpp

clean:
	rm -rf *o main