CC = mpic++
CFLAGS = -O3 --std=c++11 

all: mpi_main.exe

mpi_main.exe: main.o random.o city.h path.h salesman.h
	$(CC) random.o main.o -o mpi_main.exe $(CFLAGS)

main.o: main.cpp city.h path.h salesman.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o *.exe

copy:
	scp -r matteo.razza@tolab.fisica.unimi.it:/home/studenti/matteo.razza/NSL_Esercitazione_10/path* /Users/TeoRazza/Library/Mobile\ Documents/com~apple~CloudDocs/Università\ degli\ Studi\ di\ Milano/Triennale\ \(Fisica\ L30\)/III\ Anno/II\ Semestre/Simulazione\ Numerica/Esercitazione_10/RESULTS/08_CORE_INDIPENDENT