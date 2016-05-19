#include <iostream>
#include <string>
#include <vector>
#include "Proteine.hpp"
#include "AcideAmine.hpp"
#include "gensvg.hpp"

//using namespace std

int main() {

//Example 3
   std::string s = "HHPHPPPPHHPHHPPPHPHHPPH";

//Example 2
//   std::string s = "HHHPPHPHPHPPHPHPHPPH";

//Example 1
//   std::string s = "PHPPHPHPPPPHPPPPPHPPHHP";
   Proteine protein1 = Proteine(s);
   
   std::cout << "La longeur de la proteine est : " << protein1.l << std::endl;
   
   std::cout << "La sequence de la proteine est : ";
   
   for(int i = 0; i < protein1.l; i++){
         std::cout << protein1.proteine[i]->valeur;
   }
               
   std::cout << std::endl;

   std::cout << "Les acides animes polaires sont aux positions : ";
   
   for(int i = 0; i < protein1.l; i++){
         std::cout << protein1.proteine[i]->indice << " ";
   }
   
   std::cout << std::endl;
   
   protein1.Ranger();

   std::cout << "neff : " << protein1.neff << std::endl;
   
   showProtein(protein1);
   
   for(int i = 0; i < protein1.l; i++){
         std::cout << "Indice : " << protein1.proteine[i]->indice << " ";
         std::cout << protein1.proteine[i]->valeur << " ";
         std::cout << "Position : (" << protein1.proteine[i]->x << "," << protein1.proteine[i]->y << ") ";
         std::cout << std::endl;
   }
     
   
}

 
