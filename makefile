COMP=mpiCC

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
FLAGS=-lstdc++ -O3 -Werror
else
FLAGS=-O3 -Werror
endif

EXEC=./prj/darius
SRC_DIR=./src/
SRCS:=$(shell find $(SRC_DIR)*.cpp)
OBJ_DIR= ./obj/
OBJS:= $(SRCS:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)
LIB_DIR= ./lib/
LIB= $(LIB_DIR)libMithra.a

all: $(EXEC)

full: clean all

lib: all $(LIB)

$(OBJ_DIR)%.o: $(SRC_DIR)%.cpp
	@mkdir -p $(@D)
	$(COMP) $(FLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	$(COMP) $(FLAGS) $^ -o $@

$(LIB): $(OBJS)
	@mkdir -p $(@D)
	ar rcs $@ $^

clean:
	rm -rf $(EXEC) $(OBJ_DIR) $(LIB_DIR)
