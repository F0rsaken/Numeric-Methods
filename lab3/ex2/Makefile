MAIN_FILE=main
LIB_PATH=../../lib/
LIB_NAME=matrix-lib
LIB_PATH_2=../../lib/
LIB_NAME_2=equations-lib

all: main

main: $(MAIN_FILE).o $(LIB_NAME).o $(LIB_NAME_2).o
	g++ $(MAIN_FILE).o $(LIB_NAME).o $(LIB_NAME_2).o -o main

$(MAIN_FILE).o: $(MAIN_FILE).cpp
	g++ -c $(MAIN_FILE).cpp

$(LIB_NAME).o: $(LIB_PATH)$(LIB_NAME).cpp
	g++ -c $(LIB_PATH)$(LIB_NAME).cpp

$(LIB_NAME_2).o: $(LIB_PATH_2)$(LIB_NAME_2).cpp
	g++ -c $(LIB_PATH_2)$(LIB_NAME_2).cpp

clean:
	rm -rf *o main