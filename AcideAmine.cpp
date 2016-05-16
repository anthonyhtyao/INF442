#include <iostream>
#include "AcideAmine.hpp"

AcideAmine::AcideAmine(int i, char c, int xx, int yy) {
   
   indice = i;
   valeur = c;
   x = xx;
   y = yy;
   
}

AcideAmine::AcideAmine() {
   indice = 0;
   valeur = 'H';
   x = 0;
   y = 0;
}
