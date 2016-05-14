#ifndef GENSVGH
#define GENSVGH

#include <iostream>
#include <fstream>
#include "Proteine.hpp"
#include "AcideAnime.hpp"
#include <vector>

using namespace std;

void header(ofstream& myfile);
void footer(ofstream& myfile);
void line(int x1, int y1, int x2, int y2, ofstream& myfile);
void lien(int x1, int y1, int x2, int y2, ofstream& myfile);
void circle(int x, int y, char v, ofstream& myfile);
void showProtein(Proteine p);

#endif
