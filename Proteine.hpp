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
       
       void calculV();
       void calculVInv();
       
       void RangerAutoRight(AcideAmine* a, AcideAmine* b);
       void RangerAutoLeft(AcideAmine* a, AcideAmine* b);
       
       int nRefK();
       int calculeNeff();
       void Ranger();
       void translation();    
       //Protein sequence
       std::string sequence;
       std::vector<AcideAmine*> proteine;
       
       //Protein's length
       int l;
       int neff; 
       
       //Vector of polar Amino Acids
       std::vector<AcideAmine*> polaires;
       std::vector<AcideAmine*> typePL;
       std::vector<AcideAmine*> typePR;
       
       std::vector<std::vector<int> > v;
       std::vector<std::vector<int> > vInv;
       
};

#endif
