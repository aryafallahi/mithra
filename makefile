COMP=mpiCC

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
FLAGS=-lstdc++ -O3 -Werror
else
FLAGS=-O3 -Werror
endif

EXEC=./prj/MITHRA
SRC_DIR=./src/
SRCS:=$(shell find $(SRC_DIR)*.cpp)
HDRS:=$(shell find $(SRC_DIR)*.h)
OBJ_DIR= ./obj/
OBJS:= $(SRCS:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)
LIB_DIR= ./lib/
LIB= $(LIB_DIR)mithra.a
INC_DIR= ./include/
INC_SUBDIR= $(INC_DIR)mithra/
INCS:= $(HDRS:$(SRC_DIR)%.h=$(INC_SUBDIR)%.h)

all: $(EXEC)

full: clean all

install: all $(LIB) $(INCS)

$(OBJ_DIR)%.o: $(SRC_DIR)%.cpp
	@mkdir -p $(@D)
	$(COMP) $(FLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	$(COMP) $(FLAGS) $^ -o $@

$(LIB): $(OBJS)
	@mkdir -p $(@D)
	ar rcs $@ $^

$(INC_SUBDIR)%.h: $(SRC_DIR)%.h
	@mkdir -p $(@D)
	ln -s $< $@

clean:
	rm -rf $(EXEC) $(OBJ_DIR) $(LIB_DIR) $(INC_DIR)
