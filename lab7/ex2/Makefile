MAIN_FILE=main
LIB_PATH=../../lib/
LIB_NAME=matrix-lib
LIB_PATH_2=../../lib/
LIB_NAME_2=diff-equations-lib
LIB_PATH_3=../../lib/
LIB_NAME_3=interpolation-lib

all: main

main: $(MAIN_FILE).o $(LIB_NAME).o $(LIB_NAME_2).o $(LIB_NAME_3).o
	g++ $(MAIN_FILE).o $(LIB_NAME).o $(LIB_NAME_2).o $(LIB_NAME_3).o -o main

$(MAIN_FILE).o: $(MAIN_FILE).cpp
	g++ -c $(MAIN_FILE).cpp

$(LIB_NAME).o: $(LIB_PATH)$(LIB_NAME).cpp
	g++ -c $(LIB_PATH)$(LIB_NAME).cpp

$(LIB_NAME_2).o: $(LIB_PATH_2)$(LIB_NAME_2).cpp
	g++ -c $(LIB_PATH_2)$(LIB_NAME_2).cpp

$(LIB_NAME_3).o: $(LIB_PATH_3)$(LIB_NAME_3).cpp
	g++ -c $(LIB_PATH_3)$(LIB_NAME_3).cpp

clean:
	rm -rf *o main