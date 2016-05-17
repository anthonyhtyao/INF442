#include <iostream>
#include <string>
#include <vector>
#include "Proteine.hpp"
#include "AcideAmine.hpp"

//using namespace std;

Proteine::Proteine(std::string s) {

    sequence = s;
    l = s.size();
    
    for(int i=0; i<l; i++){
       AcideAmine* a = new AcideAmine(i,s[i],0,i);
       proteine.push_back(a);
       if(s[i] == 'P') {
           polaires.push_back(a);
        }
    }
    
    neff = 0;
    
    calculV();
    calculVInv();


/*    for(std::string::iterator it = s.begin(); it!=s.end(); it++) {
       AcideAmine a = AcideAmine()
    }
*/                
}
 
void Proteine::calculV() {
      
      int courantI = 0;
      int courantP = 0;
      
      for(int k=0; k<l; k++){
         std::vector<int> tmp;
         if(proteine[k]->valeur == 'P') {
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
      for(  int k=0; k<l; k++){
         vInv.push_back(aux);
      }
      
      for(int k=l-1; k>=0; k--){
         if(proteine[k]->valeur == 'P') {
            if(k%2 != 0) courantI += 1;
            else courantP +=1;
         }         
         vInv[k][0] = courantI;
         vInv[k][1] = courantP;
      }
}

 
void Proteine::RangerAutoRight(AcideAmine* a, AcideAmine* b) {
   
   int i = a->indice;
   int j = b->indice;
   
   int y1 = a->y;
   int y2 = b->y; 
   int dist_y = y1 - y2 -1;
   
   int nb = j - i -1;
   
   int dist_x = (nb - dist_y)/2;
   
   for(  int k=i+1; k<=i+dist_x; k++){
      proteine[k]->x = a->x + k - i;
      proteine[k]->y = a->y;
   }
   
   for(  int k=i+dist_x+1; k<=i+2*dist_x; k++){
      proteine[k]->x = a->x + dist_x - (k - i -dist_x -1);
      proteine[k]->y = a->y -1;
   }
   
   for(  int k=j-dist_y; k<=j-1; k++){
      proteine[k]->x = a->x;
      proteine[k]->y = b->y + (j-k);
   }
   
}

int Proteine::nRefK() {
   
   int indMin1 = 0;
   int indMax1 = 0;
   
   int indMin2 = 0;
   int indMax2 = 0;
   
   int courantMax1 = 0;
   int courantMax2 = 0;
   
   for(  int k=0; k<l; k++){
      
      int min1 = std::min(v[k][0],vInv[k][1]);
      int min2 = std::min(v[k][1],vInv[k][0]);
      
      if(min1 >= courantMax1) {
         indMax1 = k;
         if(min1 > courantMax1) {
            indMin1 = k;
            courantMax1 = min1;
         }
      }
      
      if(min2 >= courantMax2) {
         indMax2 = k;
         if(min2 > courantMax2) {
            indMin2 = k;
            courantMax2 = min2;
         }
      }
   }
   
   int indMin = indMin1;
   int indMax = indMax1;
   
   if(courantMax2 > courantMax1) {
      indMin = indMin2;
      indMax = indMax2;
   }
   
//   std::cout << indMin << " " << indMax << std::endl;
   
   return (indMin + indMax)/2;
   
}

void Proteine::RangerAutoLeft(AcideAmine* a, AcideAmine* b) {
   
   int i = a->indice;
   int j = b->indice;
   
   int y1 = a->y;
   int y2 = b->y; 
   int dist_y = y1 - y2 -1;
   
   int nb = i - j -1;
   
   int dist_x = (nb - dist_y)/2;
   
   for(  int k=i-1; k >= i-dist_x; k--){
      proteine[k]->x = a->x - (i - k);
      proteine[k]->y = a->y;
   }
   
   for(  int k= i - dist_x - 1; k >= i - 2*dist_x; k--){
      
      proteine[k]->x = a->x - (2*dist_x + 1 - (i-k));
      proteine[k]->y = a->y - 1;
   }

   for(  int k=j+1; k<=j+dist_y; k++){
      proteine[k]->x = a->x;
      proteine[k]->y = b->y + k -j;
   }
  
}

void Proteine::Ranger() {
   
   bool isImpair = false;
   
   int k = nRefK();
   std::cout << "La valeur de k : " << k << std::endl;
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
   for(unsigned int r=0; r<polaires.size(); r++) {
         int w = polaires[r]->indice;
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
   for(unsigned int r=0; r<polaires.size(); r++) {
         int w = polaires[r]->indice;
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
   
            
   proteine[k]->x = 0;
   proteine[k]->y = 0;
   
   proteine[k+1]->x = 1;
   proteine[k+1]->y = 0;
   
   int tt = k - left[0];
   proteine[left[0]]->x = 0;
   proteine[left[0]]->y = -tt;
   
   proteine[right[0]]->x = 1;
   proteine[right[0]]->y = -tt;
   
   RangerAutoLeft(proteine[k],proteine[left[0]]);
   RangerAutoRight(proteine[k+1],proteine[right[0]]);
   
   for(int r = 1; r < nref; r++){
      
      tt = std::min(left[r-1] - left[r], right[r] - right[r-1]);
      proteine[left[r]]->x = 0;
      proteine[left[r]]->y = proteine[left[r-1]]->y - tt;
      
      proteine[right[r]]->x = 1;
      proteine[right[r]]->y = proteine[right[r-1]]->y - tt;
      
      RangerAutoLeft(proteine[left[r-1]], proteine[left[r]]);
      RangerAutoRight(proteine[right[r-1]],proteine[right[r]]);
      
   }
   
   proteine[0]->x = proteine[k]->x;
   proteine[0]->y = proteine[left[nref -1]]->y - left[nref-1];
   
   RangerAutoLeft(proteine[left[nref-1]],proteine[0]);   
   
   proteine[l-1]->x = proteine[k+1]->x;
   proteine[l-1]->y = proteine[right[nref -1]]->y - (l-1 - right[nref -1]);   
   
   RangerAutoRight(proteine[right[nref-1]], proteine[l-1]);

   neff = calculeNeff();
   
   translation();
/*
   for(int r=0; r<nref; r++){
      int ind = leftaux[nref-1-r];
      typePL.push_back(proteine[ind]);
      ind = right[r];
      typePR.push_back(proteine[ind]);
   }
*/   
}

void Proteine::translation() {
   for (int i = 0; i < l; i++) {
      proteine[i]->x += l;
      proteine[i]->y = proteine[i]->y*(-1) + l;
   }
}
   
int Proteine::calculeNeff() {
   int n = 0;
   for (unsigned  int i = 0; i < polaires.size(); i++) {
      int ind = polaires[i]->indice;
      AcideAmine* a1 = new AcideAmine();
      a1 = polaires[i];
      for (unsigned  int j = i+1; j < polaires.size(); j++) {
         int ind1 = polaires[j]->indice;
         if (ind + 1 != ind1) {
            AcideAmine* a2 = new AcideAmine();
            a2 = proteine[ind1];
            if (a1->x == a2->x) {
               if (a1->y == a2->y + 1 || a1->y == a2->y - 1) {
                  n++;
                  typePL.push_back(a1);
                  typePR.push_back(a2);
               }
            }
            else if (a1->y == a2->y) {
               if (a1->x == a2->x + 1 || a1->x == a2->x - 1) {
                  n++;
                  typePL.push_back(a1);
                  typePR.push_back(a2);
               }
            }
         }
      }
   }
   return n;
}

bool Proteine::notOverlap(int i) {
   bool b = true;
   for (int k = 0; k <i; k++) {
      b = b && ((proteine[k]->x != proteine[i]->x) || (proteine[k]->y != proteine[i]->y));   
   }
   return b;
}

void Proteine::RangerRecursif(int i, Proteine p) {
   if(i == l) {}
   else {
      int xPre = p.proteine[i-1]->x;
      int yPre = p.proteine[i-1]->y;
      
      p.proteine[i]->x = xPre;
      p.proteine[i]->y = yPre + 1;
      
      if(p.notOverlap(i)) {
         RangerRecursif(i+1, p);
         
         if(p.test()) {
            p.neff = p.calculeNeff();
            if(p.neff > neff) {
               for(int q = 0; q<l ; q++){
                  *proteine[q] = *p.proteine[q];
               }
            neff = p.neff;
            }
         }
      }
       
      p.proteine[i]->x = xPre;
      p.proteine[i]->y = yPre - 1;
      
      if(p.notOverlap(i)) {
        
         RangerRecursif(i+1, p);
         
         if(p.test()) {
            p.neff = p.calculeNeff();
            if(p.neff > neff) {
               for(int q = 0; q<l ; q++){
                  *proteine[q] = *p.proteine[q];
               }
            neff = p.neff;
            }
         }
      }
      
      p.proteine[i]->x = xPre - 1;
      p.proteine[i]->y = yPre;
      
      if(p.notOverlap(i)) {
         
         RangerRecursif(i+1, p);
        
         if(p.test()) {
            p.neff = p.calculeNeff(); 
            if(p.neff > neff) {
               for(int q = 0; q<l ; q++){
                  *proteine[q] = *p.proteine[q];
               }            
            neff = p.neff;
            }
         }
      }
      
      p.proteine[i]->x = xPre + 1;
      p.proteine[i]->y = yPre;
      
      if(p.notOverlap(i)) {
         
         RangerRecursif(i+1, p);
         
         if(p.test()) {
            p.neff = p.calculeNeff();
            if(p.neff > neff) {
               for(int q = 0; q<l ; q++){
                  *proteine[q] = *p.proteine[q];
               }
            neff = p.neff;
           }
         }        
      }
   }
}

bool Proteine::test() {
   
   for (int i=0; i<l; i++) {
      for (int j = i+1; j<l;j++) {
         if ((proteine[i]->x == proteine[j]->x)&&(proteine[i]->y == proteine[j]->y)) return false;   
      }   
   }   
   return true;
}
