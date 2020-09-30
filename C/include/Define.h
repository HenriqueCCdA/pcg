#ifndef _DEFINE_
  #define _DEFINE_
/*...*/
  #include<stdlib.h>
  #include<stdio.h>
/*...................................................................*/

/*... solver*/
  #define PCG        1
/*...................................................................*/

/*... macro para acesso matricial em vetores*/
  #define   MAT2D(i,j,vector,col)           vector[i*col+j]
  #define   MAT3D(i,j,k,vector,col1,col2)   vector[i*col1*col2+col2*j+k]
/*...................................................................*/

/*... definicao de funcoes*/
  #define min(a, b)  (((a) < (b)) ? (a) : (b))
  #define max(a, b)  (((a) > (b)) ? (a) : (b))
  #define vectorPlusOne(v,n,i)  for(i=0;i<n;i++) v[i]++ 
  #define vectorMinusOne(v,n,i) for(i=0;i<n;i++) v[i]--  
/*...................................................................*/

/*...*/
  #define HccaAbs(a)  (((a) < 0.0e0) ? (-a) : (a))
/*...................................................................*/


/*... Saida de Erro*/                                                  
#define ERRO_RCM fprintf(stderr,"\nrcm - fatal error!\n")

  #define ERRO_OP(line,file,func,op)\
    fprintf(stderr,"Opecao %d e invalida!!\n",op);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
          ,file,func,line);\
    exit(EXIT_FAILURE);
  
  #define ERRO_NORM(line,file,func,lin,col)\
    fprintf(stderr,"Erro no calulo da norma!!\n");\
    fprintf(stderr,"Numero de linha  : %d\n"\
                   "Numero de colunas: %d\n",lin,col);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
          ,file,func,line);\
    exit(EXIT_FAILURE);
 
  #define ERRO_GERAL(file,func,str)\
    fprintf(stderr,"Erro: %s!!\n",str);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\n",file,func);\
    exit(EXIT_FAILURE);

  #define ERRO_MALLOC(point,str,line,file,func)\
     if(point == NULL){\
     fprintf(stderr,"Erro na alocacao do vetor %s\n",str);\
     fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
            ,file,func,line);\
     exit(EXIT_FAILURE);}

/*...................................................................*/
#endif/*_DEFINE_*/
