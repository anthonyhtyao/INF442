#ifndef PROTEINEH
#define PROTEINEH

#include <iostream>
#include <string>
#include <vector>

//using namespace std;

class Proteine {
   
    public:
    
       //Constructor
       Proteine(std::string s);
    
       //Protein sequence
       std::string sequence;
       std::vector<char> p;
       
       //Protein's length
       int l; 
       
       //Protein's structure
       //For each Amino Acid, we stock its position compared with its previous
       //Up, down, left and right are noted respectively by 1,2,3,4 
       std::vector<int> position;
       
};

#endif
