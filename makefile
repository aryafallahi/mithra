SHELL = /bin/sh
COMP = mpic++

CFLAGS+=-std=c++11

EXEC = prj/MITHRA
SRC_DIR = src
SRCS := $(shell find $(SRC_DIR)/*.cpp)
HDRS := $(shell find $(SRC_DIR)/*.h)
OBJ_DIR = obj
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
LIB_DIR = lib
LIB = $(LIB_DIR)/libmithra.a
INC_DIR = include/mithra/

.PHONY: all clean debug install

all: $(EXEC)

debug: CFLAGS += -g
debug: all

install: all $(LIB) install-incs

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

clean:
	$(RM) -r $(EXEC) $(OBJ_DIR) $(LIB_DIR) $(INC_DIR)
