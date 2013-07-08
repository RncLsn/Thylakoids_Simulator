CC= gcc
LFLAG= -L/usr/local/lib -lgsl -lgslcblas -lm
CFLAG= -O2  -Wall -I./head #-DDEBUG_ROOT #-g  #-DDEBUG_DEBYE#-DDEBUG_EQ#-DDEBUG_ROOT#-DDEBUG # -DDEBUG -DDEBUG_ROOT   -DDEBUG -DDEBUG_1
HEAD= ./head
SRC= ./src
BIN= ./bin
RES= ./results

all: main

main: $(BIN)/main.o $(BIN)/io_util.o $(BIN)/strutture_dati.o $(BIN)/calcolo_util.o 
	$(CC) $(BIN)/main.o $(LFLAG) $(BIN)/io_util.o $(BIN)/strutture_dati.o $(BIN)/calcolo_util.o -o main

$(BIN)/main.o: $(SRC)/main.c $(HEAD)/io_util.h $(HEAD)/strutture_dati.h  $(HEAD)/calcolo_util.h
	$(CC) -c $(SRC)/main.c -o $(BIN)/main.o $(CFLAG)

$(BIN)/io_util.o: $(SRC)/io_util.c $(HEAD)/io_util.h
	$(CC) -c $(SRC)/io_util.c -o $(BIN)/io_util.o $(CFLAG)

$(BIN)/strutture_dati.o: $(SRC)/strutture_dati.c $(HEAD)/strutture_dati.h
	$(CC) -c $(SRC)/strutture_dati.c -o $(BIN)/strutture_dati.o $(CFLAG)

$(BIN)/calcolo_util.o: $(SRC)/calcolo_util.c $(HEAD)/calcolo_util.h
	$(CC) -c $(SRC)/calcolo_util.c -o $(BIN)/calcolo_util.o $(CFLAG)

# run + plot
.phony: run
run:
	./main
	gnuplot "gnuplot_script"
	xdg-open $(RES)/input.png
	xdg-open $(RES)/time-length.png 


#clean binaries and results
.phony: clean_all
clean_all:
	rm $(BIN)/*.o main $(RES)/*

#clean only results
.phony: clean_results
clean_results:
	rm $(RES)/*

#debug
.phony: memcheck
memcheck:
	valgrind ./main

# only charts
.phony: plot_do
plot_do:
	gnuplot "gnuplot_script"

