CC = g++
CFLAGS = -O3 --std=c++11 

all: main.exe

main.exe: main.o random.o city.h path.h salesman.h
	$(CC) random.o main.o -o main.exe $(CFLAGS)

main.o: main.cpp city.h path.h salesman.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm *.o *.exe
