#ifndef _ERROR_H

/*...*/
  #include<stdlib.h>
  #include<stdio.h>
/*...................................................................*/

/*... Saida de Erro*/                                                  
  #define ERRO_MALLOC(point,str,line,file,func)\
     if(point == NULL){\
     fprintf(stderr,"Erro na alocacao do vetor %s\n",str);\
     fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
            ,file,func,line);\
     exit(EXIT_FAILURE);}

  #define ERRO_FILE(point,str,line,file,func)\
     if(point == NULL){\
     fprintf(stderr,"Erro na abertura do arquivo %s\n",str);\
     fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
            ,file,func,line);\
     exit(EXIT_FAILURE);}

/*...................................................................*/
#endif