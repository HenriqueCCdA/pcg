#ifndef _DEFINE_
  #define _DEFINE_
/*...*/
  #include<stdlib.h>
  #include<stdio.h>
/*...................................................................*/

/*...*/
  #ifdef _MSC_VER
/*  #define LONG_INT __int64*/
    typedef __int64 LONG_INT;
    #define RESTRICT __restrict
  #else
    #define LONG_INT long  
/*  typedef __int64 LONG_INT;*/
    #define RESTRICT restrict
  #endif
/*...................................................................*/

/*...*/
  #define NONE        0
  #define DIAGONAL    1
  #define DILU        2
  #define ILU         3  
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

#endif/*_DEFINE_*/
