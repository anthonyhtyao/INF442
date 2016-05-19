#include <iostream>
#include <string>
#include <vector>
#include "Proteine.hpp"
#include "AcideAmine.hpp"
#include "gensvg.hpp"

int main() {

   //std::string s = "HPPHPPHPPHPH";
   std::string s = "PPPPHHPHHPPHPPPHPP";
   //std::string s = "HPHPPHHPHPPHPHHPPHPH";
   //std::string s = "HHPPPPHPPHHPPHHHHH";
   //std::string s = "HHPHPPPPHHPHHPPPHPHHPPH";
   Proteine* protein = new Proteine(s);
   Proteine* p = new Proteine(s);
   std::vector<int> end;
   for (int i=0; i<p->l; i++){
      if (i ==1) end.push_back(1);
      else end.push_back(0);
   }
   std::cout << "La longeur de la proteine est : " << protein->l << std::endl;
   
   
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
   
   //protein->RangerRecursif(1,p);
   std::cout << "Solution approche, to find neff seuil" << std::endl;
   protein->Ranger();
   std::cout << "Solution approche end" << std::endl;
   int nbOpt = protein->RangerAll(p,end);
   
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
   std::cout << "Le nombre de solutions optimales : " << nbOpt << std::endl;
   protein->calculeNeff();
   protein->translation();
   showProtein(*protein);
   
}
