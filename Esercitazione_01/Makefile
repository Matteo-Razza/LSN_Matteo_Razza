CC = g++
CFLAGS = -Wall -O3 --std=c++11

# Specifica i nomi degli eseguibili e i relativi file sorgente
EXECUTABLES = main1 main2 main3 main4 main5 main6
SOURCES = $(EXECUTABLES:=.cpp)

all: $(EXECUTABLES)

# Regola generica per la compilazione degli eseguibili
%: %.cpp random.o
	$(CC) -o $@ $^ $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm -f $(EXECUTABLES) *.o seed.out