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
         std::string sCourant = s.substr(a*(r-1),r-1);
         int c = (int) sCourant[r-1-1]-48;
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
      sendbuffer[t] = (int) s[t] - 48;
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


   //std::string sPro = "HPPHPPHPPHPH";
   //std::string sPro = "PPPPHHPHHPPHPPPHPP";
   std::string sPro = "HPHPPHHPHPPHPHHPPHPH";
   Proteine* protein = new Proteine(sPro);
   Proteine* p = new Proteine(sPro);
   std::vector<int> end;
   for(int i=0; i<p->l;i++) end.push_back(0);
   //end[1] = 1;
   int ind;
   int nbOpt;
   if(taskid < puiss){
      for(int i=rang-2; i>=0; i--){
         if(rbuf[i] != 3){
            ind = i;
            break;
         }
      }
      for(int i = 0; i<(rang-1); i++){
         /*if(i<ind) end[i] = rbuf[i];
         else if (i == ind) end[i] = rbuf[i]+1;*/
         p->pos[i+2] = rbuf[i];
         if(rbuf[i] == 0){
            p->proteine[i+2]->x = p->proteine[i+1]->x +1;
            p->proteine[i+2]->y = p->proteine[i+1]->y;
            if(i == rang-2) end[i+2] = rbuf[i+2]+1;
         }
         else if(rbuf[i] == 1){
            p->proteine[i+2]->x = p->proteine[i+1]->x;
            p->proteine[i+2]->y = p->proteine[i+1]->y +1;
            if(i == rang-2) end[i+2] = rbuf[i+2]+1;
         }
         else if(rbuf[i] == 2){
            p->proteine[i+2]->x = p->proteine[i+1]->x -1;
            p->proteine[i+2]->y = p->proteine[i+1]->y;
            if(i == rang-2) end[i+2] = rbuf[i+2]+1;
         }
         else{
            p->proteine[i+2]->x = p->proteine[i+1]->x;
            p->proteine[i+2]->y = p->proteine[i+1]->y-1;
            if(i == rang-2) end[i+1] = 1;
         }
      }
      for (int i = rang+1; i<p->l;i++) {
         p->proteine[i]->x = p->proteine[i-1]->x +1;
         p->proteine[i]->y = p->proteine[i-1]->y;
      }
      //protein->RangerRecursif(rang,p);
      protein->Ranger();
      nbOpt = protein->RangerAll(p,end);
   }

   int neffLocal = protein->neff;
   int neffMax;
   MPI_Reduce(&neffLocal, &neffMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
   int nbOptLocal = 0;
   //if(neffLocal == neffMax) nbOptLocal = protein->nbOpt;
   if(neffLocal == neffMax) nbOptLocal = nbOpt;
   int nbOptGlobal;
   MPI_Reduce(&nbOptLocal, &nbOptGlobal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);   
   //std::cout << "La valeur de neff locale est : " << neffLocal << std::endl;
   //std::cout << "Le nombre de solutions optimales est : " << nbOptLocal << std::endl;
   //nbOptGlobal = nbOptGlobal /2;
   nbOptGlobal = nbOptGlobal;
   if(taskid == 0){
      std::cout << "La valeur de neff maximale est : " << neffMax << std::endl;
      std::cout << "Le nombre de solutions optimales est : " << nbOptGlobal << std::endl;
   }
   MPI_Finalize();
   return 0;
}
