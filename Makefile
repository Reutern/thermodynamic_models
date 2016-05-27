CC = g++ -std=c++11 -fopenmp -O3 -mtune=native -march=native #-g #-O0 #for gdb, valgrind etc

GSL_DIR = usr/lib
EIGEN_DIR = /usr/include/eigen3/  
CMAES_DIR = /opt/libcmaes/src

ODIR=.obj

INC = -I$(GSL_DIR) -I$(EIGEN_DIR) -I$(CMAES_DIR)
CFLAGS = -O3 $(INC) 

LIB_DIRS = -L$(GSL_DIR)/lib -L$(EIGEN_DIR) -L$(CMAES_DIR)/lib
LIBS = -lgsl -lgslcblas -lm -lcmaes
LFLAGS = -O3 $(LIB_DIRS) $(LIBS)

BIN = seq2expr

all: $(BIN)

clean:
	rm -f $(BIN)
	rm -f $(ODIR)/*.o

_OBJ = Tools.o SeqAnnotator.o ExprPredictor.o OccPredictor.o type.o seq2expr.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

DEPS = Tools.h SeqAnnotator.h ExprPredictor.h OccPredictor.h type.h param.h

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

seq2expr : $(OBJ)
	$(CC) -o $@ $(OBJ) $(LFLAGS)

.PHONY: clean


