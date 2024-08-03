# Compiler
CC = gcc
NVCC = nvcc

# Flags
CFLAGS = -c
NVCCFLAGS = -arch=sm_75

# Directories
LIB_DIR = src/lib
SRC_DIR = src

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

# Object files for different optimization levels
OBJS_SEQ_O0 = $(LIB_DIR)/stack_O0.o $(LIB_DIR)/queue_O0.o $(LIB_DIR)/graph_O0.o $(LIB_DIR)/state_O0.o $(SRC_DIR)/$(EXEC_SEQ_BASE)_O0.o
OBJS_SEQ_O1 = $(LIB_DIR)/stack_O1.o $(LIB_DIR)/queue_O1.o $(LIB_DIR)/graph_O1.o $(LIB_DIR)/state_O1.o $(SRC_DIR)/$(EXEC_SEQ_BASE)_O1.o
OBJS_SEQ_O2 = $(LIB_DIR)/stack_O2.o $(LIB_DIR)/queue_O2.o $(LIB_DIR)/graph_O2.o $(LIB_DIR)/state_O2.o $(SRC_DIR)/$(EXEC_SEQ_BASE)_O2.o
OBJS_SEQ_O3 = $(LIB_DIR)/stack_O3.o $(LIB_DIR)/queue_O3.o $(LIB_DIR)/graph_O3.o $(LIB_DIR)/state_O3.o $(SRC_DIR)/$(EXEC_SEQ_BASE)_O3.o

OBJS_PAR_O0 = $(LIB_DIR)/stack_O0.o $(LIB_DIR)/graph_O0.o $(LIB_DIR)/state_O0.o
OBJS_PAR_O1 = $(LIB_DIR)/stack_O1.o $(LIB_DIR)/graph_O1.o $(LIB_DIR)/state_O1.o
OBJS_PAR_O2 = $(LIB_DIR)/stack_O2.o $(LIB_DIR)/graph_O2.o $(LIB_DIR)/state_O2.o
OBJS_PAR_O3 = $(LIB_DIR)/stack_O3.o $(LIB_DIR)/graph_O3.o $(LIB_DIR)/state_O3.o

# Rule for compiling files
all: $(EXEC_SEQ_O0) $(EXEC_SEQ_O1) $(EXEC_SEQ_O2) $(EXEC_SEQ_O3) $(EXEC_PAR_O0) $(EXEC_PAR_O1) $(EXEC_PAR_O2) $(EXEC_PAR_O3)

# Command for object files using gcc with optimization levels
$(LIB_DIR)/%_O0.o: $(LIB_DIR)/%.c
	$(CC) $(CFLAGS) -O0 $< -o $@

$(LIB_DIR)/%_O1.o: $(LIB_DIR)/%.c
	$(CC) $(CFLAGS) -O1 $< -o $@

$(LIB_DIR)/%_O2.o: $(LIB_DIR)/%.c
	$(CC) $(CFLAGS) -O2 $< -o $@

$(LIB_DIR)/%_O3.o: $(LIB_DIR)/%.c
	$(CC) $(CFLAGS) -O3 $< -o $@

# Command for object file of exec_seq using gcc with optimization levels
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

# Cleaning
clean:
	rm -f $(LIB_DIR)/*_O0.o $(LIB_DIR)/*_O1.o $(LIB_DIR)/*_O2.o $(LIB_DIR)/*_O3.o \
	$(SRC_DIR)/*_O0.o $(SRC_DIR)/*_O1.o $(SRC_DIR)/*_O2.o $(SRC_DIR)/*_O3.o \
	$(EXEC_SEQ_O0) $(EXEC_SEQ_O1) $(EXEC_SEQ_O2) $(EXEC_SEQ_O3) \
	$(EXEC_PAR_O0) $(EXEC_PAR_O1) $(EXEC_PAR_O2) $(EXEC_PAR_O3)
