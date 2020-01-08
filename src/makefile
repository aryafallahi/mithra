COMP=mpiCC

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
FLAGS=-lstdc++ -O3 -Werror
else
FLAGS=-O3 -Werror
endif

EXEC=../darius
SRC_DIR=./
SRCS:=$(shell find $(SRC_DIR)*.cpp)
OBJ_DIR= ../obj/
OBJS:= $(SRCS:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)

all: $(EXEC)	

full: clean all

$(OBJ_DIR)%.o: $(SRC_DIR)%.cpp
	@mkdir -p $(@D)
	$(COMP) $(FLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	$(COMP) $(FLAGS) $^ -o $@

clean:
	rm -rf $(EXEC) $(OBJ_DIR)