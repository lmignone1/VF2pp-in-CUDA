# Compiler
CC = gcc
NVCC = nvcc

# Flags
CFLAGS = -c
NVCCFLAGS = -arch=sm_75

# Directories
LIB_DIR = src/lib
SRC_DIR = src
TEST_DIR = test

# Target
EXEC_SEQ_BASE = vf2pp_sequential
EXEC_PAR_BASE = vf2pp_parallel

# Executables with optimization levels
EXEC_SEQ_O0 = $(EXEC_SEQ_BASE)_O0
EXEC_SEQ_O1 = $(EXEC_SEQ_BASE)_O1
EXEC_SEQ_O2 = $(EXEC_SEQ_BASE)_O2
EXEC_SEQ_O3 = $(EXEC_SEQ_BASE)_O3

EXEC_PAR_O0 = $(EXEC_PAR_BASE)_O0
EXEC_PAR_O1 = $(EXEC_PAR_BASE)_O1
EXEC_PAR_O2 = $(EXEC_PAR_BASE)_O2
EXEC_PAR_O3 = $(EXEC_PAR_BASE)_O3

EXEC_SEQ_TEST1 = $(EXEC_SEQ_BASE)_test1
EXEC_PAR_TEST1 = $(EXEC_PAR_BASE)_test1

EXEC_SEQ_TEST2 = $(EXEC_SEQ_BASE)_test2
EXEC_PAR_TEST2 = $(EXEC_PAR_BASE)_test2

EXEC_SEQ_TEST3 = $(EXEC_SEQ_BASE)_test3
EXEC_PAR_TEST3 = $(EXEC_PAR_BASE)_test3

# Object files for different optimization levels
OBJS_SEQ_O0 = $(LIB_DIR)/stack_O0.o $(LIB_DIR)/queue_O0.o $(LIB_DIR)/graph_O0.o $(LIB_DIR)/state_O0.o $(SRC_DIR)/$(EXEC_SEQ_BASE)_O0.o
OBJS_SEQ_O1 = $(LIB_DIR)/stack_O1.o $(LIB_DIR)/queue_O1.o $(LIB_DIR)/graph_O1.o $(LIB_DIR)/state_O1.o $(SRC_DIR)/$(EXEC_SEQ_BASE)_O1.o
OBJS_SEQ_O2 = $(LIB_DIR)/stack_O2.o $(LIB_DIR)/queue_O2.o $(LIB_DIR)/graph_O2.o $(LIB_DIR)/state_O2.o $(SRC_DIR)/$(EXEC_SEQ_BASE)_O2.o
OBJS_SEQ_O3 = $(LIB_DIR)/stack_O3.o $(LIB_DIR)/queue_O3.o $(LIB_DIR)/graph_O3.o $(LIB_DIR)/state_O3.o $(SRC_DIR)/$(EXEC_SEQ_BASE)_O3.o

OBJS_PAR_O0 = $(LIB_DIR)/stack_O0.o $(LIB_DIR)/graph_O0.o $(LIB_DIR)/state_O0.o
OBJS_PAR_O1 = $(LIB_DIR)/stack_O1.o $(LIB_DIR)/graph_O1.o $(LIB_DIR)/state_O1.o
OBJS_PAR_O2 = $(LIB_DIR)/stack_O2.o $(LIB_DIR)/graph_O2.o $(LIB_DIR)/state_O2.o
OBJS_PAR_O3 = $(LIB_DIR)/stack_O3.o $(LIB_DIR)/graph_O3.o $(LIB_DIR)/state_O3.o

OBJS_SEQ_TEST1 = $(LIB_DIR)/stack.o $(LIB_DIR)/queue.o $(LIB_DIR)/graph.o $(LIB_DIR)/state.o $(TEST_DIR)/$(EXEC_SEQ_TEST1).o
OBJS_SEQ_TEST2 = $(LIB_DIR)/stack.o $(LIB_DIR)/queue.o $(LIB_DIR)/graph.o $(LIB_DIR)/state.o $(TEST_DIR)/$(EXEC_SEQ_TEST2).o
OBJS_SEQ_TEST3 = $(LIB_DIR)/stack.o $(LIB_DIR)/queue.o $(LIB_DIR)/graph.o $(LIB_DIR)/state.o $(TEST_DIR)/$(EXEC_SEQ_TEST3).o

OBJS_PAR_TEST = $(LIB_DIR)/stack.o $(LIB_DIR)/graph.o $(LIB_DIR)/state.o

# Rule for compiling files
all: $(EXEC_SEQ_O0) $(EXEC_SEQ_O1) $(EXEC_SEQ_O2) $(EXEC_SEQ_O3) $(EXEC_PAR_O0) $(EXEC_PAR_O1) $(EXEC_PAR_O2) $(EXEC_PAR_O3)

# Command for object files using gcc with optimization levels. We are compiling the library files with different optimization levels
$(LIB_DIR)/%_O0.o: $(LIB_DIR)/%.c
	$(CC) $(CFLAGS) -O0 $< -o $@

$(LIB_DIR)/%_O1.o: $(LIB_DIR)/%.c
	$(CC) $(CFLAGS) -O1 $< -o $@

$(LIB_DIR)/%_O2.o: $(LIB_DIR)/%.c
	$(CC) $(CFLAGS) -O2 $< -o $@

$(LIB_DIR)/%_O3.o: $(LIB_DIR)/%.c
	$(CC) $(CFLAGS) -O3 $< -o $@

# Command for object file of sequential version using gcc with optimization levels
$(SRC_DIR)/$(EXEC_SEQ_BASE)_O0.o: $(SRC_DIR)/$(EXEC_SEQ_BASE).c
	$(CC) $(CFLAGS) -O0 $< -o $@

$(SRC_DIR)/$(EXEC_SEQ_BASE)_O1.o: $(SRC_DIR)/$(EXEC_SEQ_BASE).c
	$(CC) $(CFLAGS) -O1 $< -o $@

$(SRC_DIR)/$(EXEC_SEQ_BASE)_O2.o: $(SRC_DIR)/$(EXEC_SEQ_BASE).c
	$(CC) $(CFLAGS) -O2 $< -o $@

$(SRC_DIR)/$(EXEC_SEQ_BASE)_O3.o: $(SRC_DIR)/$(EXEC_SEQ_BASE).c
	$(CC) $(CFLAGS) -O3 $< -o $@

# Commands for sequential executable with different optimizations
$(EXEC_SEQ_O0): $(OBJS_SEQ_O0)
	$(CC) $(OBJS_SEQ_O0) -o $@

$(EXEC_SEQ_O1): $(OBJS_SEQ_O1)
	$(CC) $(OBJS_SEQ_O1) -o $@

$(EXEC_SEQ_O2): $(OBJS_SEQ_O2)
	$(CC) $(OBJS_SEQ_O2) -o $@

$(EXEC_SEQ_O3): $(OBJS_SEQ_O3)
	$(CC) $(OBJS_SEQ_O3) -o $@

# Commands for parallel executable with different optimizations
$(EXEC_PAR_O0): $(OBJS_PAR_O0)
	$(NVCC) $(NVCCFLAGS) -O0 $(SRC_DIR)/$(EXEC_PAR_BASE).cu $(OBJS_PAR_O0) -o $@

$(EXEC_PAR_O1): $(OBJS_PAR_O1)
	$(NVCC) $(NVCCFLAGS) -O1 $(SRC_DIR)/$(EXEC_PAR_BASE).cu $(OBJS_PAR_O1) -o $@

$(EXEC_PAR_O2): $(OBJS_PAR_O2)
	$(NVCC) $(NVCCFLAGS) -O2 $(SRC_DIR)/$(EXEC_PAR_BASE).cu $(OBJS_PAR_O2) -o $@

$(EXEC_PAR_O3): $(OBJS_PAR_O3)
	$(NVCC) $(NVCCFLAGS) -O3 $(SRC_DIR)/$(EXEC_PAR_BASE).cu $(OBJS_PAR_O3) -o $@

# Test rule
test: $(EXEC_SEQ_TEST1) $(EXEC_PAR_TEST1) $(EXEC_SEQ_TEST2) $(EXEC_PAR_TEST2) $(EXEC_SEQ_TEST3) $(EXEC_PAR_TEST3)

# Command for object files using gcc. We are compiling the library files
$(LIB_DIR)/%.o: $(LIB_DIR)/%.c
	$(CC) $(CFLAGS) $< -o $@

# Command for object test file of sequential version using gcc
$(TEST_DIR)/$(EXEC_SEQ_TEST1).o: $(TEST_DIR)/$(EXEC_SEQ_TEST1).c
	$(CC) $(CFLAGS) $< -o $@

$(TEST_DIR)/$(EXEC_SEQ_TEST2).o: $(TEST_DIR)/$(EXEC_SEQ_TEST2).c
	$(CC) $(CFLAGS) $< -o $@

$(TEST_DIR)/$(EXEC_SEQ_TEST3).o: $(TEST_DIR)/$(EXEC_SEQ_TEST3).c
	$(CC) $(CFLAGS) $< -o $@

# Command for executable test of sequential version
$(EXEC_SEQ_TEST1): $(OBJS_SEQ_TEST1)
	$(CC) $(OBJS_SEQ_TEST1) -o $@

$(EXEC_SEQ_TEST2): $(OBJS_SEQ_TEST2)
	$(CC) $(OBJS_SEQ_TEST2) -o $@

$(EXEC_SEQ_TEST3): $(OBJS_SEQ_TEST3)
	$(CC) $(OBJS_SEQ_TEST3) -o $@

# Command for executable test of parallel version
$(EXEC_PAR_TEST1): $(OBJS_PAR_TEST)
	$(NVCC) $(NVCCFLAGS) $(TEST_DIR)/$(EXEC_PAR_TEST1).cu $(OBJS_PAR_TEST) -o $@

$(EXEC_PAR_TEST2): $(OBJS_PAR_TEST)
	$(NVCC) $(NVCCFLAGS) $(TEST_DIR)/$(EXEC_PAR_TEST2).cu $(OBJS_PAR_TEST) -o $@

$(EXEC_PAR_TEST3): $(OBJS_PAR_TEST)
	$(NVCC) $(NVCCFLAGS) $(TEST_DIR)/$(EXEC_PAR_TEST3).cu $(OBJS_PAR_TEST) -o $@

# Cleaning
clean:
	rm -f $(LIB_DIR)/*.o $(SRC_DIR)/*.o $(TEST_DIR)/*.o \
	$(EXEC_SEQ_O0) $(EXEC_SEQ_O1) $(EXEC_SEQ_O2) $(EXEC_SEQ_O3) \
	$(EXEC_PAR_O0) $(EXEC_PAR_O1) $(EXEC_PAR_O2) $(EXEC_PAR_O3) \
	$(EXEC_SEQ_TEST1) $(EXEC_PAR_TEST1) $(EXEC_SEQ_TEST2) $(EXEC_PAR_TEST2) \
	$(EXEC_SEQ_TEST3) $(EXEC_PAR_TEST3)
