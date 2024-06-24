CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo
all: simulator.exe
simulator.exe : system.o NSL_SIMULATOR.o particle.o random.o
	$(CC) system.o particle.o random.o NSL_SIMULATOR.o -o simulator.exe $(AFLAGS)
NSL_SIMULATOR.o : NSL_SIMULATOR.cpp system.h particle.h
	$(CC) -c NSL_SIMULATOR.cpp -o NSL_SIMULATOR.o $(CFLAGS)
system.o : system.cpp system.h random.h particle.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)
particle.o : particle.cpp particle.h
	$(CC) -c particle.cpp -o particle.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
liquid:
	mv ../OUTPUT/potential_energy.dat ../RESULTS/Esercitazione_07/liquid_potential_energy_NVE.dat
	mv ../OUTPUT/gofr.dat ../RESULTS/Esercitazione_07/liquid_gofr_NVE.dat
	mv ../OUTPUT/pressure.dat ../RESULTS/Esercitazione_07/liquid_pressure_NVE.dat
copygas:
	scp razza@159.149.2.20:/home/razza/NSL_SIMULATOR/SSH/gas* /Users/TeoRazza/Library/Mobile\ Documents/com~apple~CloudDocs/Università\ degli\ Studi\ di\ Milano/Triennale\ \(Fisica\ L30\)/III\ Anno/II\ Semestre/Simulazione\ Numerica/Esercitazione_04-06-07/NSL_SIMULATOR/RESULTS/Esercitazione_07