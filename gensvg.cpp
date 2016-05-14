#include <iostream>
#include <fstream>
#include "Proteine.hpp"
#include "AcideAmine.hpp"
#include <vector>

using namespace std;

// Constant declarations specifying the dimensions of the image and
// the rectangular bar.

// Prolog for the SVG image

void header(ofstream& myfile, int IMG_WIDTH, int IMG_HEIGHT) {
  myfile << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl
         << "<svg width=\"" << IMG_WIDTH 
         << "\" height=\"" << IMG_HEIGHT << "\">" << endl
         << "<title>Proteine</title>" << endl;
}

// Epilog for the SVG image

void footer(ofstream& myfile) {
  myfile << "</svg>" << endl;
}

// Generation a line
void line(int x1, int y1, int x2, int y2, ofstream& myfile) {
    myfile << "<line x1=\"" << x1 <<"\" y1=\"" << y1 <<"\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" style=\"stroke:black;stroke-width:1\" />" << endl;
  
}

// Generation a dotted line
void lien(int x1, int y1, int x2, int y2, ofstream& myfile) {
    myfile << "<line stroke-dasharray=\"2 1\" x1=\"" << x1 <<"\" y1=\"" << y1 <<"\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" style=\"stroke:black;stroke-width:1\" />" << endl;
  
}
// Generation a circle

void circle(int x, int y, char v, ofstream& myfile) {
    if (v == 'P')
      myfile << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"6\" fill=\"black\" />" << endl;
    else
      myfile << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"6\" stroke=\"black\" stroke-width=\"1\" fill=\"white\" />" << endl;
}

// Read input and generate SVG image

void showProtein(Proteine p) {
  int IMG_WIDTH = p.l*40;
  int IMG_HEIGHT = p.l*40;
  int coeff = 20;
  ofstream myfile;
  myfile.open ("example.svg");
  header(myfile, IMG_WIDTH, IMG_HEIGHT);
  for (unsigned int i=0; i< p.typePL.size(); i++) {
    cout << p.typePL[i].y << endl;
    lien(p.typePL[i].x*coeff, p.typePL[i].y*coeff, p.typePR[i].x*coeff, p.typePR[i].y*coeff, myfile);
  }
  for (unsigned int i = 0; i < p.l; i++) {
    vector<AcideAmine> lst = p.proteine;
    if (i != p.l-1) {
      line(lst[i].x*coeff, lst[i].y*coeff, lst[i+1].x*coeff, lst[i+1].y*coeff, myfile);
    }
    circle(lst[i].x*coeff, lst[i].y*coeff, lst[i].valeur, myfile);
  }
  footer(myfile);
  myfile.close();
}
