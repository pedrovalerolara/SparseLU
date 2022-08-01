include make.inc

SRC		 =	../src
BIN		 =	BIN
INC		 +=	-I../include
MMSPARSE =  ../mmsparse

all:

	#DDSS routines	

	#TILE TEST
	$(CC) $(INC) $(FLAGS) -c $(SRC)/ddss_tile.c -o $(BIN)/ddss_tile.o $(LIBS)
	$(CC) $(INC) $(FLAGS) -c $(SRC)/ddss_flat2tiled.c -o \
          $(BIN)/ddss_flat2tiled.o $(LIBS)
	$(CC) $(INC) $(FLAGS) -c $(SRC)/ddss_tiled2flat.c -o \
          $(BIN)/ddss_tiled2flat.o $(LIBS)
	$(CC) $(INC) $(FLAGS) $(BIN)/ddss_tile.o $(BIN)/ddss_flat2tiled.o  \
      	  $(BIN)/ddss_tiled2flat.o ddss_test_dtile.c -o test_dtile $(LIBS)

	#BLAS3 routines
	
	#DGEMM TEST
	$(CC) $(INC) $(FLAGS) -c $(SRC)/ddss_dgemm.c -o \
	      $(BIN)/ddss_dgemm.o $(LIBS)  
	$(CC) $(INC) $(FLAGS) -c $(SRC)/kdgemm.c -o $(BIN)/kdgemm.o $(LIBS) 
	$(CC) $(INC) $(FLAGS) $(BIN)/ddss_tile.o \
          $(BIN)/ddss_flat2tiled.o $(BIN)/ddss_tiled2flat.o \
		  $(BIN)/ddss_dgemm.o $(BIN)/kdgemm.o \
		  ddss_test_dgemm.c -o test_dgemm $(LIBS)
	
	#LAPACK routines

	#DNPGETRF TEST
	@echo beginning dnpgetrf test compile

	@echo compiling dsss_dnpgetrf.c
	$(CC) $(INC) $(FLAGS) -c $(SRC)/dsss_dnpgetrf.c -o \
          $(BIN)/dsss_dnpgetrf.o $(LIBS)
	
	@echo compiling dnpgetrf.c
	$(CC) $(INC) $(FLAGS) -c $(SRC)/dnpgetrf.c -o \
         $(BIN)/dnpgetrf.o $(LIBS)

	@echo compiling kdnpgetrf.c
	$(CC) $(INC) $(FLAGS) -c $(SRC)/kdnpgetrf.c -o \
		$(BIN)/kdnpgetrf.o $(LIBS)

	@echo compiling dsss_flat2tiled.c
	$(CC) $(INC) $(FLAGS) -c $(SRC)/dsss_flat2tiled.c -o \
		$(BIN)/dsss_flat2tiled.o $(LIBS)

	@echo compiling dsss_tiled2flat.c
	$(CC) $(INC) $(FLAGS) -c $(SRC)/dsss_tiled2flat.c -o \
		$(BIN)/dsss_tiled2flat.o $(LIBS)

	$(CC) $(INC) $(FLAGS) $(BIN)/ddss_tile.o \
		$(BIN)/kdnpgetrf.o \
	  	$(BIN)/ddss_flat2tiled.o \
		$(BIN)/ddss_tiled2flat.o \
	  	$(BIN)/dsss_flat2tiled.o \
		$(BIN)/dsss_tiled2flat.o \
           	$(BIN)/ddss_dgemm.o \
		$(BIN)/kdgemm.o \
  		$(BIN)/dsss_dnpgetrf.o \
		$(BIN)/dnpgetrf.o \
            dsss_test_dnpgetrf.c -o test_sparse_dnpgetrf $(LIBS)	


clean:
	rm BIN/*
	rm test_* 
