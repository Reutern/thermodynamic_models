CC = g++ -std=c++11 -fopenmp -O3 -mtune=native -march=native #-g #-O0 #for gdb, valgrind etc

GSL_DIR = usr/local

INC = -I$(GSL_DIR)
CFLAGS = -O3 $(INC) 

LIB_DIRS = -L$(GSL_DIR)/lib
LIBS = -lgsl -lgslcblas -lm
LFLAGS = -O3 $(LIB_DIRS) $(LIBS)

BIN = seq2expr 

all: $(BIN)

clean:
	rm -f $(BIN)
	rm -f *.o

Tools.o : Tools.h Tools.cpp
	$(CC) $(CFLAGS) -c Tools.cpp
siman.o : siman.h siman.cpp
	$(CC) $(CFLAGS) -c siman.cpp
SeqAnnotator.o : Tools.h SeqAnnotator.h param.h SeqAnnotator.cpp
	$(CC) $(CFLAGS) -c SeqAnnotator.cpp
ExprPredictor.o : Tools.h siman.h SeqAnnotator.h ExprPredictor.h param.h ExprPredictor.cpp 
	$(CC) $(CFLAGS) -c ExprPredictor.cpp
OccPredictor.o : OccPredictor.h param.h SeqAnnotator.h ExprPredictor.h OccPredictor.cpp 
	$(CC) $(CFLAGS) -c OccPredictor.cpp
seq2expr.o : Tools.h siman.h SeqAnnotator.h ExprPredictor.h OccPredictor.h param.h seq2expr.cpp
	$(CC) $(CFLAGS) -c seq2expr.cpp

seq2expr : Tools.o siman.o SeqAnnotator.o ExprPredictor.o OccPredictor.o seq2expr.o 
	$(CC) -o $@ Tools.o siman.o SeqAnnotator.o ExprPredictor.o OccPredictor.o seq2expr.o $(LFLAGS)

