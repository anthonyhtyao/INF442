#include <iostream>
#include <fstream>
#include "Proteine.hpp"
#include "AcideAnime.hpp"
#include <vector>

using namespace std;

// Constant declarations specifying the dimensions of the image and
// the rectangular bar.

const unsigned IMG_WIDTH = 70;
const unsigned IMG_HEIGHT = 80;
const unsigned IMG_MARGIN = 10;

// Prolog for the SVG image

void header(ofstream& myfile) {
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
    myfile << "<line x1=\"" << x1 <<"\" y1=\"" << y1 <<"\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" style=\"stroke:rgb(255,0,0);stroke-width:2\" />" << endl;
  
}

// Generation a dotted line
void lien(int x1, int y1, int x2, int y2, ofstream& myfile) {
    myfile << "<line stroke-dasharray=\"1 1\" x1=\"" << x1 <<"\" y1=\"" << y1 <<"\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" style=\"stroke:rgb(255,0,0);stroke-width:2\" />" << endl;
  
}
// Generation a circle

void circle(int x, int y, char v, ofstream& myfile) {
    if (v == 'P')
      myfile << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"1.5\" fill=\"black\" />" << endl;
    else
      myfile << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"1.5\" fill=\"white\" />" << endl;
}

// Read input and generate SVG image

void showProtein(Proteine p) {
  ofstream myfile;
  myfile.open ("example.svg");
  header(myfile);
  for (unsigned int i = 0; i < p.l; i++) {
    vector<AcideAnime> lst = p.proteine;
    circle(lst[i].x*5 + 50, -lst[i].y*5 + 10, lst[i].valeur, myfile);
    if (i != p.l-1) {
      line(lst[i].x*5 + 50, -lst[i].y*5 + 10, lst[i+1].x*5 + 50, -lst[i+1].y*5 + 10, myfile);
    }
  }
  for (unsigned int i=0; i< p.typePL.size(); i++) {
    cout << p.typePL[i].y << endl;
    lien(p.typePL[i].x*5 + 50, -(p.typePL[i].y*5) + 10, p.typePR[i].x*5 + 50, -(p.typePR[i].y*5) + 10, myfile);
  }
  footer(myfile);
  myfile.close();
/*
  unsigned m, n;

  cin >> m;
  cin >> n;

  header();
  line();
  bar(m, n);
  footer();
*/
}
