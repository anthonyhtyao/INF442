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
       // For ranger()
       void RangerAutoRight(AcideAmine* a, AcideAmine* b);
       void RangerAutoLeft(AcideAmine* a, AcideAmine* b);
       // Return index k which has nref max
       int nRefK();
       // Find solution approche
       void Ranger();
       void translation();

       // Return true if no acide anime overlaps before i
       bool notOverlap(int i);
       // Return true if no up before i
       bool firstUp(int i);
       int nextPosition(int ind, bool r);
       bool test();
       // Calcule neff for this protein
       int calculeNeff();

       // Find solution optimal by using recursif function
       int RangerRecursif(int i, Proteine* p);

       // Find next structure possible
       bool shift(int seuil);
       // Find solution optimal by using loop
       int RangerAll(Proteine* p, std::vector<int> end); 
          
       //Protein sequence
       std::string sequence;
       std::vector<AcideAmine*> proteine;

       //Store relative position of i by i-1. 0,1,2,3 present right,up,left,down respectively
       std::vector<int> pos;
       std::vector<int> contactPossible; 
       //Protein's length
       int l;
       int neff; 
       int nbOpt; 
       
       //Vector of polar Amino Acids
       std::vector<AcideAmine*> hydrophobes;
       std::vector<AcideAmine*> typeHL;
       std::vector<AcideAmine*> typeHR;
       
       std::vector<std::vector<int> > v;
       std::vector<std::vector<int> > vInv;
       
};

#endif
