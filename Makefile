# Compilatore
CC = gcc

# Flags
CFLAGS = -c
LDFLAGS =

# Directory dei sorgenti e degli oggetti
SRC_DIR = src/lib
SEQ_DIR = src/sequential

# File sorgenti
SRC_FILES = stack.c queue.c graph.c state.c
SEQ_FILES = main.c

# File oggetto
OBJS = $(SRC_DIR)/stack.o $(SRC_DIR)/queue.o $(SRC_DIR)/graph.o $(SRC_DIR)/state.o $(SEQ_DIR)/main.o

# Target eseguibile
EXEC = main

# Regola per compilare tutti i file oggetto
all: $(EXEC)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $< -o $@

$(SEQ_DIR)/main.o: $(SEQ_DIR)/main.c
	$(CC) $(CFLAGS) $< -o $@

# Regola per creare l'eseguibile
$(EXEC): $(OBJS)
	$(CC) $(OBJS) -o $@

# Pulizia
clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all clean
