#include "AcideAmine.hpp"
#include "Proteine.hpp"
#include <iostream>
#include "mpi.h"
#include <string>
#include <stdlib.h>

std::string genere(int r, int puiss){
   if(r==1) {
      return "013";
   }
   else{
      std::string result = "";
      std::string s = genere(r-1, puiss/3);
      for(int a=0; a<puiss/3; a++){
         std::string sCourant = s.substr(a*(r-1),(a+1)*(r-1)-1);
         int c = (int) sCourant[r-1-1];
         if(c == 0){
            result = result + sCourant + "0";
            result = result + sCourant + "1";
            result = result + sCourant + "3";
         }
         if(c == 1){
            result = result + sCourant + "0";
            result = result + sCourant + "1";
            result = result + sCourant + "2";
         }
         if(c == 2){
            result = result + sCourant + "1";
            result = result + sCourant + "2";
            result = result + sCourant + "3";
         }
         if(c == 3){
            result = result + sCourant + "0";
            result = result + sCourant + "2";
            result = result + sCourant + "3";
         }
      }
      return result;
   }
}

int main(int argc, char **argv) {
   const int root = 0;

   int numtasks, taskid;
   MPI_Status status;

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
   MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

   int rang = 1;
   int puiss = 1;
   while(puiss <= numtasks){
      puiss = puiss * 3;
      rang += 1;
   }

   rang -= 1;
   puiss = puiss / 3;

   int *sendbuffer;
   int *displs, *scounts, rbuf[rang-1];
   sendbuffer = (int *) malloc((rang-1)*numtasks*(sizeof(int)));
   std::string s = genere(rang-1,puiss);
   for(unsigned int t=0; t<s.size(); t++){
      sendbuffer[t] = (int) s[t];
   }
   for(unsigned int t = s.size(); t < numtasks*(rang-1); t++){
      sendbuffer[t] = 0;
   }

   displs = (int *)malloc(numtasks*(sizeof(int)));
   scounts = (int *)malloc(numtasks*(sizeof(int)));
   for(int i=0; i<numtasks; ++i){
      displs[i] = i*(rang-1);
      scounts[i] = rang-1;
   }

   MPI_Scatterv(sendbuffer, scounts, displs, MPI_INT, rbuf, (rang-1), MPI_INT,
         0, MPI_COMM_WORLD);


   std::string sPro = "HPPHPPHPPHPH";
   Proteine* protein = new Proteine(sPro);
   Proteine* p = new Proteine(sPro);

   if(taskid < puiss){

      for(int i = 0; i<(rang-1); i++){
         if(rbuf[i] == 0){
            p->proteine[i+2]->x = p->proteine[i+1]->x +1;
            p->proteine[i+2]->y = p->proteine[i+1]->y;
         }
         if(rbuf[i] == 1){
            p->proteine[i+2]->x = p->proteine[i+1]->x;
            p->proteine[i+2]->y = p->proteine[i+1]->y +1;
         }
         if(rbuf[i] == 2){
            p->proteine[i+2]->x = p->proteine[i+1]->x -1;
            p->proteine[i+2]->y = p->proteine[i+1]->y;
         }
         if(rbuf[i] == 3){
            p->proteine[i+2]->x = p->proteine[i+1]->x;
            p->proteine[i+2]->y = p->proteine[i+1]->y-1;
         }
      }

      protein->RangerRecursif(rang,p);
   }

   int neffLocal = protein->neff;
   int neffMax;
   MPI_Reduce(&neffLocal, &neffMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

   if(taskid == 0){
      std::cout << "La valeur de neff maximale est : " << neffMax << std::endl;
   }

   MPI_Finalize();
   return 0;
}
