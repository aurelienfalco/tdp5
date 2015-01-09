CC = mpicc
EX = mpiexec
EXEC = lu
ALG = dger dgetf2 lu dscal util dgemm_scalaire dtrsm dgesv dgetrf
SRC=$(ALG:=.c)
OBJ=$(SRC:.c=.o)
CFLAGS = -std=c99 -g -O0 -Wall -Wextra
n = 2
m = 128
N = 5
e = 0.001
p = 0
seq = 0
stat = sp-proc

all: $(EXEC)

exec: $(EXEC)
	$(EX) -np $(n) $(EXEC) $(m) $(seq) $(p)

qsub: $(EXEC)
	rm -rf res.*
	@qsub batch; 
	@sleep 3; 
	@cat res.*

$(EXEC): $(OBJ) main.c
	$(CC) $(CFLAGS) $^ -o $@ -lm

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -lm

tests: $(OBJ) tests.o
	$(CC) $(CFLAGS) $^ -o $@ -lm

test: tests
	./tests $(N) $(e) $(p)

stat: $(EXEC)
	./stats.sh

plot: 
	@gnuplot -e "xname='Nombre de processus';name='$(stat)';data='$(stat).data';output='$(stat).png" plot_sp.gp 
	@eog $(stat).png 2>/dev/null &

plot-sp: 
	@gnuplot -e "xname='Longueur du côté des matrices';name='Speedup';data='sp-size.data';output='sp-size.png" plot_sp.gp 
	@eog sp-size.png 2>/dev/null &

clean:
	rm -rf *.o $(EXEC) tests *~

