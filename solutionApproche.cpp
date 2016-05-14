#include <iostream>
#include <string>
#include <vector>
#include "Proteine.hpp"
#include "AcideAmine.hpp"
#include "gensvg.hpp"

//using namespace std

int main() {
   
   std::string s = "PPHPHHHHPPHPPHHHPHPPHHP";
   Proteine protein1 = Proteine(s);
   
   std::cout << "La longeur de la proteine est : " << protein1.l << std::endl;
   
   std::cout << "La sequence de la proteine est : ";
   
   for(std::vector<AcideAmine>::iterator it = protein1.proteine.begin();
         it != protein1.proteine.end(); it++){
         std::cout << it->valeur;
   }
               
   std::cout << std::endl;
   
   std::cout << "Les acides animes polaires sont aux positions : ";
   
   for(std::vector<AcideAmine>::iterator it = protein1.polaires.begin();
         it != protein1.polaires.end(); it++){
         std::cout << it->indice << " ";
   }
   
   std::cout << std::endl;
   
   protein1.Ranger();

   showProtein(protein1);
   
   for(std::vector<AcideAmine>::iterator it = protein1.proteine.begin();
         it != protein1.proteine.end(); it++){
         std::cout << "Indice : " << it->indice << " ";
         std::cout << it->valeur << " ";
         std::cout << "Position : (" << it->x << "," << it->y << ") ";
         std::cout << std::endl;
   }
     
   
}

 
