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
  double dot(double *RESTRICT x1,double *RESTRICT x2,int const n);
/*...................................................................*/

/*level 2*/
/* ... matriz cheia*/
  void matVecFull( double *RESTRICT a
                 , double *RESTRICT x
                 , double *RESTRICT y
                 , int const nLin , int const nCol
                 , short const code);
/*...................................................................*/
#endif
