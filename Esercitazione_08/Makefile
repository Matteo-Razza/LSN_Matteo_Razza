CC = g++
CFLAGS = -O3 --std=c++11

one: main.exe

two: sa.exe

sa.exe: system.o sa.o random.o
	$(CC) system.o random.o sa.o -o sa.exe $(CFLAGS)

sa.o: sa.cpp system.h
	$(CC) -c sa.cpp -o sa.o $(CFLAGS)

main.exe: system.o main.o random.o
	$(CC) system.o random.o main.o -o main.exe $(CFLAGS)

main.o: main.cpp system.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

system.o: system.cpp system.h random.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm -f *.o main.exe sa.exe

remove:
	rm -f *.dat
