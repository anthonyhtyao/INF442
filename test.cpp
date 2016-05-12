#include <iostream>
#include <string>
#include <vector>
#include "Proteine.hpp"

//using namespace std

int main() {
   
   std::string s = "HPPPHPPHPH";
   Proteine protein1 = Proteine(s);
   
   for(std::vector<char>::iterator it = protein1.p.begin();
         it != protein1.p.end(); it++){
         std::cout << *it;
   }
   
   std::cout << std::endl;
   
}

