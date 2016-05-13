CC = g++

Proteine.o: Proteine.cpp Proteine.hpp
	$(CC) -c Proteine.cpp

test.o: test.cpp Proteine.hpp
	$(CC) -c test.cpp

test: test.o Proteine.o
	$(CC) Proteine.o test.o -o test

clean:
	rm -f test
	rm -f *.o
