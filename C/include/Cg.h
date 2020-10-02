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
  #include"HccaStdBool.h"
  #include"Error.h"
/*...................................................................*/
 
/*========================= Iterativos ==============================*/
/*... gradiente conjugado precondicionado*/
  void pcg(int const neq       , double *RESTRICT a
          , double *RESTRICT b , double *RESTRICT x 
          , double const tol   , char preC
          , unsigned int maxit , bool newX          
          , FILE* fileSolvLog  , bool log
          , void(*matvec)()    , double(*dot)());
/*...................................................................*/

#endif/*_SOLV_H*/
