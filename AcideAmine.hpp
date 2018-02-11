#ifndef ACIDEAMINEH
#define ACIDEAMINEH

#include <iostream>

class AcideAmine {

   public:
      
      //Constructor
      AcideAmine(int indice, char valeur, int x, int y);
      AcideAmine();
      ~AcideAmine() {};
      
      //Indice de l'acide anime
      int indice;
      
      //Abscisse et ordonnee
      int x;
      int y;
      
      //H ou P
      char valeur;
};
      
#endif      

