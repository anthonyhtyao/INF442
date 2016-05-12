#include <iostream>
#include <string>
#include <vector>
#include "Proteine.hpp"

//using namespace std;

Proteine::Proteine(std::string s) {

    sequence = s;
    l = s.size();
    
/*    for(unsigned int i=0; i<l; i++){
       proteine.push_back(s[i]);
    }
*/

    for(std::string::iterator it = s.begin(); it!=s.end(); it++) {
       p.push_back(*it);
    }

 }
