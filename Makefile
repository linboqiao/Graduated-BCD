CC = g++

CFLAGS = -O -g -Wall -fmessage-length=0 -ansi -D_GNU_SOURCE -fexceptions -DMX_COMPAT_32 -DNDEBUG -fopenmp 
CLIBS = -O -g -Wall -fmessage-length=0 -L./ -lstdc++ -fopenmp 


INC = -I ../eigen-3.2.10

SRC_DIR = .
OBJ_DIR = .
BIN_DIR = .

APP = main
SRC = $(wildcard $(SRC_DIR)/*.cpp) 
OBJ = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(notdir $(SRC)))

BIN = $(BIN_DIR)/$(APP)

all: $(BIN)
$(BIN): $(OBJ)
	$(CC) -o $@ $(OBJ) $(CLIBS)
$(OBJ_DIR)/%.o:$(SRC_DIR)/%.cpp
	$(CC) -c $< -o $@ $(INC) $(CFLAGS)


clean:
	rm  $(OBJ_DIR)/*.o 
	rm  $(BIN_DIR)/$(APP)

