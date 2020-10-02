#ifndef _PRECOND_H
  #include<stdlib.h>

  #include"Define.h"

  double* preMake(double *a, int const nEq, char preC);

  void preCondSolver(double *m     , double *a
                   , double *r    , double *x
                   , int const nEq, char preC);


  void dilu(double *m    , double *a, int const nEq);
  void fb_dilu(double *r, double *a, double *d, double *x, int const nEq);

#endif