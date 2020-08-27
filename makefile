SHELL = /bin/sh
COMP = mpic++

CFLAGS+=-std=c++11
CFLAGS+=-O3

PREFIX ?= .
LIB_DIR = $(PREFIX)/lib
INC_DIR = $(PREFIX)/include/mithra/
 
EXEC = prj/MITHRA
SRC_DIR = src
SRCS := $(shell find $(SRC_DIR)/*.cpp)
HDRS := $(shell find $(SRC_DIR)/*.h)
OBJ_DIR = obj
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
LIB = $(OBJ_DIR)/libmithra.a

.PHONY: all clean debug install

all: $(EXEC) $(LIB)

debug: CFLAGS += -g
debug: all

install: all install-incs install-lib

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(COMP) $(CFLAGS) -c $< -o $@

$(EXEC): $(OBJS)
	$(COMP) $(LDFLAGS) $^ -o $@

$(LIB): $(OBJS)
	@mkdir -p $(@D)
	$(AR) rcs $@ $^

install-incs:
	@mkdir    -p $(INC_DIR)
	install $(HDRS) $(INC_DIR)

install-lib:
	@mkdir -p $(LIB_DIR)
	install $(LIB) $(LIB_DIR)

clean:
	$(RM) $(EXEC) $(OBJS) $(LIB)
	rmdir $(OBJ_DIR)
