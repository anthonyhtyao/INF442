#include <iostream>
#include <string>
#include <vector>
#include "Proteine.hpp"
#include "AcideAmine.hpp"
#include "gensvg.hpp"

int main() {
  
   std::string s = "HHHPPHPHPHPPHPHPHPPH";
//   std::string s = "HPPHPPHPPHPHP";
//   std::string s = "PPPPHHPHHPPHPPPHPP";
//   std::string s = "PHPPHPHPPPPHPPPPPH";
   Proteine* protein = new Proteine(s);
   Proteine* p = new Proteine(s);
   
   std::cout << "La longeur de la proteine est : " << protein->l << std::endl;
   
   std::cout << "La sequence de la proteine est : ";
   
   for(int i = 0; i < protein->l; i++){
         std::cout << protein->proteine[i]->valeur;
   }
               
   std::cout << std::endl;
   std::cout << "Les acides animes hydrophobes sont aux positions : ";
   
   
   for(unsigned int i = 0; i < protein->hydrophobes.size(); i++){
         std::cout << protein->hydrophobes[i]->indice << " ";
   }
   
   std::cout << std::endl;
   std::cout << std::endl;
   std::cout << "La structure avant recherche :" << std::endl;
   
   for(int i = 0; i < protein->l; i++){
         std::cout << "Indice : " << protein->proteine[i]->indice << " ";
         std::cout << protein->proteine[i]->valeur << " ";
         std::cout << "Position : (" << protein->proteine[i]->x << "," <<
               protein->proteine[i]->y << ") ";
         std::cout << std::endl;
   }
   
   protein->RangerRecursif(2,p);
   
   
   std::cout << std::endl;
   std::cout << "La structure apres recherche : " << std::endl;
   
   for(int i = 0; i < protein->l; i++){
         std::cout << "Indice : " << protein->proteine[i]->indice << " ";
         std::cout << protein->proteine[i]->valeur << " ";
         std::cout << "Position : (" << protein->proteine[i]->x << "," <<
               protein->proteine[i]->y << ") ";
         std::cout << std::endl;
   }
   
   std::cout << "La valeur de Neff vaut : " << protein->neff << std::endl;
   protein->translation();
   showProtein(*protein);
   
}
