CC = g++ -Wall

AcideAmine.o: AcideAmine.cpp AcideAmine.hpp
	$(CC) -c AcideAmine.cpp

Proteine.o: Proteine.cpp Proteine.hpp AcideAmine.hpp
	$(CC) -c $<

gensvg.o: gensvg.cpp gensvg.hpp Proteine.hpp AcideAmine.hpp
	$(CC) -c $<

solutionApproche.o: solutionApproche.cpp Proteine.hpp gensvg.hpp AcideAmine.hpp
	$(CC) -c $<

solutionApproche: solutionApproche.o Proteine.o AcideAmine.o gensvg.o
	$(CC) -o $@ $^
	./solutionApproche
	gnome-open example.svg

clean:
	rm -f solutionApproche
	rm -f *.o
	rm -f *.svg
