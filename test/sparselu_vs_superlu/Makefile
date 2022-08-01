include make.inc

all: 
	$(CC) $(CFLAGS) -c matrix_methods.c -o matrix_methods.o
	$(CC) $(CFLAGS) main.c -o main $(SUPERLULIB) $(BLASLIB) $(BIN)/*.o matrix_methods.o
	$(CC) $(CFLAGS) matrix_gen.c -o matrix_gen matrix_methods.o

clean:
	rm main
	rm matrix_gen
