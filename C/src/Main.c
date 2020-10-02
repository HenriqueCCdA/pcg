#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include"../include/Define.h"
#include"../include/Cg.h"
#include"../include/HccaBlas.h"
#include"../include/ReadFile.h"
#include"../include/HccaStdBool.h"

int main(int argc, char *argv[]) {

  int i;
  int preC[] = {0, 1, 2};
  char *name[] = {"CG", "JCG", "DILU"};
  char nameFile[100];
  int nEq;
  double *a, *b, *x;
  FILE *fLog = NULL;
  void *matVec = &matVecFull, *inner = &dot;

/*... lendo o sistema de equacao*/
  readSystem(&a, &b, &x, &nEq);
/*...................................................................*/

/*...*/
  for(i = 0; i < 3; i++)
  {
/*... abertura do arquivo de log*/
    strcpy(nameFile, name[i]); 
    strcat(nameFile,".txt");
    
    fLog = fopen(nameFile, "w");
    ERRO_FILE(fLog, nameFile, __LINE__, __LINE__, __func__, );
/*..................................................................*/
    
/*...*/
    pcg(nEq        , a
      , b          , x 
      , 1.e-11     , i
      , 10000      , true
      , fLog       , true
      , matVec     , inner);
/*...................................................................*/
    fclose(fLog);

  }
/*...................................................................*/

/*... liberando memoria*/
  free(a);
  free(b);
  free(x);
/*...................................................................*/

  return EXIT_SUCCESS;

}