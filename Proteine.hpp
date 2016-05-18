#ifndef PROTEINEH
#define PROTEINEH

#include <iostream>
#include <string>
#include <vector>
#include "AcideAmine.hpp"

//using namespace std;

class Proteine {
   
    public:
    
       //Constructor
       Proteine(std::string s);
       ~Proteine() {};
       
       void calculV();
       void calculVInv();
       
       void RangerAutoRight(AcideAmine* a, AcideAmine* b);
       void RangerAutoLeft(AcideAmine* a, AcideAmine* b);
       
       int nRefK();
       void Ranger();
       void translation();

       // Return true if no acide anime overlaps
       bool notOverlap(int i);
       bool firstUp(int i);
       int nextPosition(int ind);
       bool test();
       int calculeNeff();
       int RangerRecursif(int i, Proteine* p);
       void shift(); 
       bool shift();
       int RangerAll(Proteine* p, std::vector<int> end); 
          
       //Protein sequence
       std::string sequence;
       std::vector<AcideAmine*> proteine;
       std::vector<int> pos; 
       //Protein's length
       int l;
       int neff; 
       
       //Vector of polar Amino Acids
       std::vector<AcideAmine*> hydrophobes;
       std::vector<AcideAmine*> typePL;
       std::vector<AcideAmine*> typePR;
       
       std::vector<std::vector<int> > v;
       std::vector<std::vector<int> > vInv;
       
};

#endif
