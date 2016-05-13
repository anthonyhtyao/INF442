#include <iostream>
#include <fstream>
#include "Proteine.hpp"
#include "AcideAnime.hpp"
#include <vector>

using namespace std;

// Constant declarations specifying the dimensions of the image and
// the rectangular bar.

const unsigned IMG_WIDTH = 100;
const unsigned IMG_HEIGHT = 100;
const unsigned IMG_MARGIN = 10;
const unsigned BAR_WIDTH = IMG_WIDTH - 2 * IMG_MARGIN;
const unsigned BAR_HEIGHT = IMG_HEIGHT - 2 * IMG_MARGIN;

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

// Generation a circle

void circle(int x, int y, char v, ofstream& myfile) {
    if (v == 'P')
      myfile << "<circle cx=\"" << x*5 << "\" cy=\"" << y*5 << "\" r=\"4\" fill=\"black\" />" << endl;
    else
      myfile << "<circle cx=\"" << x*5 << "\" cy=\"" << y*5 << "\" r=\"4\" fill=\"white\" />" << endl;
}

void bar(unsigned m, unsigned n) {
  unsigned m_width = (unsigned) (((double) m) / ((double) m + n) * BAR_WIDTH);
  unsigned n_width = BAR_WIDTH - m_width;
  cout << "<rect x=\"" << IMG_MARGIN << "\" y=\"" << IMG_MARGIN << "\""
       << " width=\"" << m_width << "\" height=\"" << BAR_HEIGHT << "\""
       << " style=\"stroke: none; fill: red;\" />" << endl;
  cout << "<rect x=\"" << IMG_MARGIN + m_width << "\" y=\"" << IMG_MARGIN << "\""
       << " width=\"" << n_width << "\" height=\"" << BAR_HEIGHT << "\""
       << " style=\"stroke: none; fill: green;\" />" << endl;
}

// Read input and generate SVG image

void showProtein(Proteine p) {
  ofstream myfile;
  myfile.open ("example.svg");
  header(myfile);
  for (vector<AcideAnime>::iterator it = p.proteine.begin(); it != p.proteine.end(); it++) {
    circle(it->x, it->y, it->valeur, myfile);
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
