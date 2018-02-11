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

solutionExhaustive.o: solutionExhaustive.cpp Proteine.hpp gensvg.hpp AcideAmine.hpp
	$(CC) -c $<

solutionExhaustiveRec: solutionExhaustive.o Proteine.o AcideAmine.o gensvg.o
	$(CC) -o $@ $^
	./solutionExhaustiveRec 1
	gnome-open example.svg

solutionExhaustiveIte: solutionExhaustive.o Proteine.o AcideAmine.o gensvg.o
	$(CC) -o $@ $^
	./solutionExhaustiveIte 0
	gnome-open example.svg

solutionDistrib: Proteine.o AcideAmine.o solutionExhaustiveDistrib.cpp
	mpic++ $^ -o $@
	mpirun -np 3 $@

clean:
	rm -f solutionApproche
	rm -f solutionExhaustive
	rm -f *.o
	rm -f *.svg
