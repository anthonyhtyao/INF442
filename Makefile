CC = g++ -Wall

AcideAnime.o: AcideAnime.cpp AcideAnime.hpp
	$(CC) -c AcideAnime.cpp

Proteine.o: Proteine.cpp Proteine.hpp AcideAnime.hpp
	$(CC) -c $<

gensvg.o: gensvg.cpp gensvg.hpp Proteine.hpp AcideAnime.hpp
	$(CC) -c $<

test.o: test.cpp Proteine.hpp gensvg.hpp AcideAnime.hpp
	$(CC) -c $<

test: test.o Proteine.o AcideAnime.o gensvg.o
	$(CC) -o $@ $^
	./test
	gnome-open example.svg

clean:
	rm -f test
	rm -f *.o
	rm -f *.svg
