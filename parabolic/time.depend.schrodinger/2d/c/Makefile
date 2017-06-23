FC=gcc
CFLAGS=-O3
LDFLAGS=-lm
EXE=tdse.x
SRC=cio.c tdse.c
OBJ=$(SRC:.c=.o)

all: $(SRC) $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ) -o $(EXE)


%.o: %.c 
	$(FC) $(CFLAGS) -c $< -o $@

clean:
	rm -f wf-*.ppm $(OBJ) $(EXE)
