#include <iostream>
#include <string>
#include <vector>
#include "Proteine.hpp"
#include "AcideAnime.hpp"

//using namespace std;

Proteine::Proteine(std::string s) {

    sequence = s;
    l = s.size();
    
    for(unsigned int i=0; i<l; i++){
       AcideAnime a = AcideAnime(i,s[i],0,i);
       proteine.push_back(a);
       if(s[i] == 'P') {
           polaires.push_back(a);
        }
    }
    
    calculV();
    calculVInv();


/*    for(std::string::iterator it = s.begin(); it!=s.end(); it++) {
       AcideAnime a = AcideAnime()
    }
*/                
}
 
void Proteine::calculV() {
      
      int courantI = 0;
      int courantP = 0;
      
      for(unsigned int k=0; k<proteine.size(); k++){
         std::vector<int> tmp;
         if(proteine[k].valeur == 'P') {
            if(k%2 != 0) courantI += 1;
            else courantP +=1;
         }
         tmp.push_back(courantI);
         tmp.push_back(courantP);
         v.push_back(tmp);
      }
}

void Proteine::calculVInv() {
      
      int courantI = 0;
      int courantP = 0;
      
      std::vector<int> aux;
      aux.push_back(0);
      aux.push_back(0);
      for(unsigned int k=0; k<l; k++){
         vInv.push_back(aux);
      }
      
      for(int k=l-1; k>=0; k--){
         if(proteine[k].valeur == 'P') {
            if(k%2 != 0) courantI += 1;
            else courantP +=1;
         }         
         vInv[k][0] = courantI;
         vInv[k][1] = courantP;
      }
}

 
void Proteine::RangerAutoRight(AcideAnime a, AcideAnime b) {
   
   int i = a.indice;
   int j = b.indice;
   
   int y1 = a.y;
   int y2 = b.y; 
   int dist_y = y1 - y2 -1;
   
   int nb = j - i -1;
   
   int dist_x = (nb - dist_y)/2;
   
   for(unsigned int k=i+1; k<=i+dist_x; k++){
      proteine[k].x = a.x + k - i;
      proteine[k].y = a.y;
   }
   
   for(unsigned int k=i+dist_x+1; k<=i+2*dist_x; k++){
      proteine[k].x = a.x + dist_x - (k - i -dist_x -1);
      proteine[k].y = a.y -1;
   }
   
   for(unsigned int k=j-dist_y; k<=j-1; k++){
      proteine[k].x = a.x;
      proteine[k].y = b.y + (j-k);
   }
   
}

int Proteine::nRefK() {
   
   int indMin = 0;
   int indMax = 0;
   int courantMax = 0;
   
   for(unsigned int k=0; k<l; k++){
      
      int min1 = std::min(v[k][0],vInv[k][1]);
      int min2 = std::min(v[k][1],vInv[k][0]);
      int max = std::max(min1,min2);
      
      if(max >= courantMax) {
         indMax = k;
         if(max > courantMax) {
            indMin = k;
            courantMax = max;
         }
      }
   }
   
//   std::cout << indMin << " " << indMax << std::endl;
   
   return (indMin + indMax)/2;
   
}

void Proteine::RangerAutoLeft(AcideAnime a, AcideAnime b) {
   
   int i = a.indice;
   int j = b.indice;
   
   int y1 = a.y;
   int y2 = b.y; 
   int dist_y = y1 - y2 -1;
   
   int nb = i - j -1;
   
   int dist_x = (nb - dist_y)/2;
   
   for(unsigned int k=i-1; k >= i-dist_x; k--){
      proteine[k].x = a.x - (i - k);
      proteine[k].y = a.y;
   }
   
   for(unsigned int k= i - dist_x - 1; k >= i - 2*dist_x; k--){
      
      proteine[k].x = a.x - (2*dist_x + 1 - (i-k));
      proteine[k].y = a.y - 1;
   }

   for(unsigned int k=j+1; k<=j+dist_y; k++){
      proteine[k].x = a.x;
      proteine[k].y = b.y + k -j;
   }
  
}

void Proteine::Ranger() {
   
   bool isImpair = false;
   
//   int k = nRefK();
   int k = 11;
   int i1 = v[k][0];
   int p1 = v[k][1];
   int i2 = vInv[k][0];
   int p2 = vInv[k][1];
  
// Get nref value 
   int nref = std::min(p1,i2);
   if(std::min(i1,p2)>std::min(p1,i2)) {
      nref = std::min(i1,p2);
      isImpair = true;
   }   
   
   std::vector<int> left;
   std::vector<int> right;
   std::vector<int> leftaux;

// Find P type Protein paired at left side  
   int compt = 0;
   std::cout << "Left P type protein ";
   for(int r=0; r<polaires.size(); r++) {
         int w = polaires[r].indice;
         if(w <= k) {
            if(isImpair && w%2 != 0){
               leftaux.push_back(w);
               compt +=1;
               std::cout << w << " ";
            }
            else if(!isImpair && w%2 ==0){
               leftaux.push_back(w);
               compt +=1;
               std::cout << w << " ";
            }
         }
         if(compt == nref) break;
   }
   std::cout << std::endl;
      
   compt = 0;

// Find P type Protein paired at right side
   std::cout << "Right P type protein ";
   for(int r=0; r<polaires.size(); r++) {
         int w = polaires[r].indice;
         if(w>k) {
            if(isImpair && w%2 == 0){
               right.push_back(w);
               compt +=1;
               std::cout << w << " ";
            }
            if(!isImpair && w%2 !=0){
               right.push_back(w);
               compt +=1;
               std::cout << w << " ";
            }
         }
         if(compt == nref) break;
   }
   std::cout << std::endl;
   
   for(int r=0; r<nref; r++){
      int ind = leftaux[nref-1-r];
      left.push_back(ind);
   }
   
/*   for(int r=polaires.size()-1; r >= 0; r--) {
         int w = polaires[r].indice;
         if(w<=k) {
            if(isImpair && w%2 != 0){
               left.push_back(w);
               compt +=1;
            }
            if(!isImpair && w%2 ==0){
               left.push_back(w);
               compt +=1;
            }
         }
         if(compt == nref) break;
  }
*/
            
   proteine[k].x = 0;
   proteine[k].y = 0;
   
   proteine[k+1].x = 1;
   proteine[k+1].y = 0;
   
   int tt = k - left[0];
   proteine[left[0]].x = 0;
   proteine[left[0]].y = -tt;
   
   proteine[right[0]].x = 1;
   proteine[right[0]].y = -tt;
   
   RangerAutoLeft(proteine[k],proteine[left[0]]);
   RangerAutoRight(proteine[k+1],proteine[right[0]]);
   
   for(int r = 1; r < nref; r++){
      
      tt = std::min(left[r-1] - left[r], right[r] - right[r-1]);
      proteine[left[r]].x = 0;
      proteine[left[r]].y = proteine[left[r-1]].y - tt;
      
      proteine[right[r]].x = 1;
      proteine[right[r]].y = proteine[right[r-1]].y - tt;
      
      RangerAutoLeft(proteine[left[r-1]], proteine[left[r]]);
      RangerAutoRight(proteine[right[r-1]],proteine[right[r]]);
      
   }
   
   proteine[0].x = proteine[k].x;
   proteine[0].y = proteine[left[nref -1]].y - left[nref-1];
   
   RangerAutoLeft(proteine[left[nref-1]],proteine[0]);   
   
   proteine[l-1].x = proteine[k+1].x;
   proteine[l-1].y = proteine[right[nref -1]].y - (l-1 - right[nref -1]);   
   
   RangerAutoRight(proteine[right[nref-1]], proteine[l-1]);
   
   for(int r=0; r<nref; r++){
      int ind = leftaux[nref-1-r];
      typePL.push_back(proteine[ind]);
      ind = right[r];
      typePR.push_back(proteine[ind]);
   }
   
}
   
      
   
   
