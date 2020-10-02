#ifndef _PRECOND_H
  #include<stdlib.h>

  #include"Define.h"
  #include"Error.h"

  double* preMake(double *RESTRICT a, int const nEq, char preC);

  void preCondSolver(double *RESTRICT m    , double *RESTRICT a
                   , double *RESTRICT r    , double *RESTRICT x
                   , int const nEq         , char preC);


  void dilu(double *RESTRICT m, double *RESTRICT a, int const nEq);
  void fb_dilu(double *RESTRICT r, double *RESTRICT a
             , double *RESTRICT d, double *RESTRICT x, int const nEq);

  void ilu(double *RESTRICT m, double *RESTRICT a, int const nEq);
  void fb_ilu(double *RESTRICT r, double *RESTRICT m
            , double *RESTRICT x, int const nEq);

#endif