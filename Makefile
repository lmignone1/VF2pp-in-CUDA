# compiler
CC = gcc

# Flags
CFLAGS = -c

# Directorries
LIB_DIR = src/lib
SRC_DIR = src

# Target
EXEC = vf2pp_sequential

# object files
OBJS = $(LIB_DIR)/stack.o $(LIB_DIR)/queue.o $(LIB_DIR)/graph.o $(LIB_DIR)/state.o $(SRC_DIR)/$(EXEC).o

# rule for compiling files
all: $(EXEC)

# command for object files
$(LIB_DIR)/%.o: $(LIB_DIR)/%.c
	$(CC) $(CFLAGS) $< -o $@

$(SRC_DIR)/$(EXEC).o: $(SRC_DIR)/$(EXEC).c
	$(CC) $(CFLAGS) $< -o $@

# command for executable
$(EXEC): $(OBJS)
	$(CC) $(OBJS) -o $@

# Cleaning
clean:
	rm -f $(OBJS) $(EXEC)