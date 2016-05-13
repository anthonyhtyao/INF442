#include <iostream>
#include <string>
#include <vector>
#include "Proteine.hpp"
#include "AcideAnime.hpp"



int main(){
   

   
   void calculVInv(proteine s) {
      
      int courantI;
      int courantP;
      
      std::vector<int> aux;
      aux.push_back(0);
      aux.push_back(0);
      for(unsigned int k=0; k<s.size(); k++){
         v.push_back(aux);
      }
      
      for(unsigned int k=s.size()-1; k>=0; k--){
         std::vector<int> tmp;
         if(s[k].valeur == 'P') {
            if(k%2 != 0) courantI += 1;
            else courantP +=1;
         }
         tmp.push_back(courantI);
         tmp.push_back(courantP);
         v[k] = tmp;
      }
   }
   

      
                     
      
