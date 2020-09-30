#ifndef _SOLV_H
  #define _SOLV_H

/*... solver*/
  #define PCG        1
/*...................................................................*/

/*...*/
  #include<stdlib.h>
  #include<math.h>
/*...................................................................*/

/*...*/
  #include"Define.h"
  #include"PreCond.h"
/*...................................................................*/
 
/*========================= Iterativos ==============================*/
/*... gradiente conjugado precondicionado*/
  void pcg(int const neq      , double *a
          , double *b          , double *x 
          , double const tol   , char preC
          , unsigned int maxit , char newX          
          , FILE* fileSolvLog  , char log
          , void(*matvec)()    , double(*dot)());
/*...................................................................*/

#endif/*_SOLV_H*/
