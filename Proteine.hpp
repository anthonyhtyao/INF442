#ifndef PROTEINEH
#define PROTEINEH

#include <iostream>
#include <string>
#include <vector>
#include "AcideAnime.hpp"

//using namespace std;

class Proteine {
   
    public:
    
       //Constructor
       Proteine(std::string s);
       
       void calculV();
       void calculVInv();
       
       void RangerAutoRight(AcideAnime a, AcideAnime b);
       void RangerAutoLeft(AcideAnime a, AcideAnime b);
       
       int nRefK();
       void Ranger();
    
       //Protein sequence
       std::string sequence;
       std::vector<AcideAnime> proteine;
       
       //Protein's length
       int l; 
       
       //Vector of polar Amino Acids
       std::vector<AcideAnime> polaires;
       
       std::vector<std::vector<int> > v;
       std::vector<std::vector<int> > vInv;
       
};

#endif
