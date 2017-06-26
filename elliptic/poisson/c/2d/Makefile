CC=gcc
CFLAGS=-O3
LDFLAGS=-lm
EXE=poisson.x
SRC=poisson.c
OBJ=$(SRC:.c=.o)

all: $(SRC) $(EXE)

$(EXE): $(OBJ)
	$(CC) $(OBJ) -o $(EXE) $(LDFLAGS)


%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXE)
