#ifndef _PRECOND_H
  #include<stdlib.h>

  #include"Define.h"

  double* preMake(double *a, int const nEq, char preC);

  void preCondSolver(double *m     , double *a
                   , double *r    , double *x
                   , int const nEq, char preC);

#endif