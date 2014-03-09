CC = g++ -fopenmp #-g for gdb

GSL_DIR = usr/include

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
SeqAnnotator.o : Tools.h SeqAnnotator.h SeqAnnotator.cpp
	$(CC) $(CFLAGS) -c SeqAnnotator.cpp
ExprPredictor.o : Tools.h SeqAnnotator.h ExprPredictor.h param.h ExprPredictor.cpp 
	$(CC) $(CFLAGS) -c ExprPredictor.cpp
seq2expr.o : Tools.h SeqAnnotator.h ExprPredictor.h seq2expr.cpp
	$(CC) $(CFLAGS) -c seq2expr.cpp

seq2expr : Tools.o SeqAnnotator.o ExprPredictor.o seq2expr.o 
	$(CC) -o $@ Tools.o SeqAnnotator.o ExprPredictor.o seq2expr.o $(LFLAGS)

