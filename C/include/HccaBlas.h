#ifndef _HCCABLAS_
  #define _HCCABLAS_
  #define HCCABLASZERO 1.0e-10
/*...*/
  #include"Define.h"
/*...*/
  #include<stdio.h>
  
/*...*/
  long flopMatVecFull(int const nLin,int const nCol);
/*...................................................................*/

/*... level 1*/
  double dot(double *x1,double *x2,int const n);
/*...................................................................*/

/*level 2*/
/* ... matriz cheia*/
  void matVecFull( double *a
                 , double *x
                 , double *y
                 , int const nLin , int const nCol
                 , short const code);
/*...................................................................*/
#endif
