#include"../include/PreCond.h"

/***************************************************************************
* DATA DE CRIACAO  : 30/09/2020                                           *
* DATA DE MODIFICAO: 01/10/2020                                           *
* ----------------------------------------------------------------------- *
* preMake: monta os precondicionadores                                    *
* ----------------------------------------------------------------------- *
* Parametros de entrada:                                                  *
* ----------------------------------------------------------------------- *
* a       - matriz de coefientes n x n                                    *
* preCond - tipo do precondicionador                                      *
*         0 - nenhum                                                      *
*         1 - Jacobi                                                      *
*         2 - DILU                                                        *
*         3 - ILU                                                         *
* neq     - numero de linhas e colunas                                    *
* ----------------------------------------------------------------------- *
* Parametros de saida:                                                    *
* ----------------------------------------------------------------------- *
* m    - retorna o inverso do coefiente da diagonal                       *
* ----------------------------------------------------------------------- *
* OBS:                                                                    *
* ----------------------------------------------------------------------- *
***************************************************************************/
double* preMake(double *RESTRICT a, int const nEq, char preC)
{
  int i;
  double *m = NULL;

/*... nenhum*/
  if(preC == NONE){
    m = (double *) malloc(sizeof(double)*nEq);
    ERRO_MALLOC(m, "m",__LINE__, __FILE__, __func__);
    for(i = 0; i < nEq; i++)
      m[i] = 1.e0;
  }
/*...................................................................*/

/*... diagonal*/
  else if(preC == DIAGONAL){
    m = (double *) malloc(sizeof(double)*nEq);
    ERRO_MALLOC(m, "m",__LINE__, __FILE__, __func__);
    for(i = 0; i < nEq; i++)
      m[i] = 1.e0/a[i*nEq + i];
  }
/*...................................................................*/

/*... DILU*/
  else if(preC == DILU){
    m = (double *) malloc(sizeof(double)*nEq);
    ERRO_MALLOC(m, "m",__LINE__, __FILE__, __func__);
    dilu(m, a, nEq);
  }
/*...................................................................*/

/*... ILU*/
  else if(preC == ILU){
    m = (double *) malloc(sizeof(double)*nEq*nEq);
    ERRO_MALLOC(m, "m",__LINE__, __FILE__, __func__);
    ilu(m, a, nEq);
  }
/*...................................................................*/

  return m;

}
/**************************************************************************/

/**************************************************************************
* DATA DE CRIACAO  : 30/09/2020                                           *
* DATA DE MODIFICAO: 01/10/2020                                           *
* ----------------------------------------------------------------------- *
* preMake: monta os precondicionadores                                    *
* ----------------------------------------------------------------------- *
* Parametros de entrada:                                                  *
* ----------------------------------------------------------------------- *
* m   -> precondicionador                                                 *
* a   -> matriz de coeficientes nxn                                       *
* r   -> matriz de coefientes n                                           *
* x   -> matriz de coefientes n                                           *
* nEq -> numero de linhas e colunas                                       *
* preC-> tipo do precondicionador                                         *
*         0 - nenhum                                                      *
*         1 - Jacobi                                                      *
*         2 - DILU                                                        *
*         3 - ILU                                                         *
* ----------------------------------------------------------------------- *
* Parametros de saida:                                                    *
* ----------------------------------------------------------------------- *
* x -> retorna a solucao                                                  *
* ----------------------------------------------------------------------- *
* OBS:                                                                    *
* ----------------------------------------------------------------------- *
***************************************************************************/
void preCondSolver(double *RESTRICT m    , double *RESTRICT a
                 , double *RESTRICT r    , double *RESTRICT x
                 , int const nEq         , char preC)
{

  int i;

/*... nenhum*/
  if(preC == NONE){
    for(i = 0; i < nEq; i++)
      x[i] = r[i];
  }

/*... diagonal*/
  else if(preC == DIAGONAL){
    for(i = 0; i < nEq; i++)
      x[i] = r[i]*m[i];
  }

/*... DILU*/
  else if(preC == DILU){
    fb_dilu(r, a, m, x, nEq);
  }

/*... ILU*/
  else if(preC == ILU){
    fb_ilu(r, m, x, nEq);
  }

}
/**************************************************************************/

/**************************************************************************
* DATA DE CRIACAO  : 30/09/2020                                           *
* DATA DE MODIFICAO: 00/00/0000                                           *
* ----------------------------------------------------------------------- *
* dilu: precondicionador DILU(0)                                          *
* ----------------------------------------------------------------------- *
* Parametros de entrada:                                                  *
* ----------------------------------------------------------------------- *
* m   -> nao definido                                                     *
* a   -> matriz de coeficientes nxn                                       *
* nEq -> numero de linhas e colunas                                       *
* ----------------------------------------------------------------------- *
* Parametros de saida:                                                    *
* ----------------------------------------------------------------------- *
* n -> precondicionador                                                   *
* ----------------------------------------------------------------------- *
* OBS:                                                                    *
* ----------------------------------------------------------------------- *
***************************************************************************/
void dilu(double *RESTRICT m, double *RESTRICT a, int const nEq)
{

  int i, j;

/*...*/
  for(i = 0; i < nEq; i++)
    m[i] = MAT2D(i, i, a, nEq);
/*........................................................................*/

/*...*/
  for(i = 0; i < nEq; i++)
    for(j = i + 1; j < nEq; j++)
      m[j] -= MAT2D(j, i, a, nEq) * MAT2D(i, j, a, nEq) / m[i];
/*........................................................................*/

}
/**************************************************************************/

/**************************************************************************
* DATA DE CRIACAO  : 30/09/2020                                           *
* DATA DE MODIFICAO: 01/10/2020                                           *
* ----------------------------------------------------------------------- *
* fb_dilu: forward e backward substituicao                                *
* ----------------------------------------------------------------------- *
* Parametros de entrada:                                                  *
* ----------------------------------------------------------------------- *
* r   - vetor de indepemdente                                             *
* a   - matriz de coeficientes                                            *
* d   - diagonal modificada do precondicionador                           *
* x   - nao definido                                                      *
* nEq - numero de equacoes                                                *
* ----------------------------------------------------------------------- *
* Parametros de saida:                                                    *
* ----------------------------------------------------------------------- *
* x    - solucao                                                          *
* ----------------------------------------------------------------------- *
* OBS:                                                                    *
* ----------------------------------------------------------------------- *
***************************************************************************/
void fb_dilu(double *RESTRICT r, double *RESTRICT a
           , double *RESTRICT d, double *RESTRICT x, int const nEq)
{

  int i, j;

/*...*/
  for(i = 0; i < nEq; i++){
    x[i] = r[i] / d[i];
    for(j = 0; j < i; j++)
      x[i] -=  MAT2D(i, j, a, nEq) * x[j] / d[i];
  }
/*........................................................................*/

/*...*/
  for(i = nEq - 1; i > -1; i--)
    for(j = i + 1; j < nEq; j++)
      x[i] -=  MAT2D(i, j, a, nEq) * x[j] / d[i];
/*........................................................................*/

}
/**************************************************************************/

/**************************************************************************
* DATA DE CRIACAO  : 02/10/2020                                           *
* DATA DE MODIFICAO: 00/00/0000                                           *
* ----------------------------------------------------------------------- *
* preMake: monta os precondicionadores                                    *
* ----------------------------------------------------------------------- *
* Parametros de entrada:                                                  *
* ----------------------------------------------------------------------- *
* m   -> nao definido                                                     *
* a   -> matriz de coeficientes nxn                                       *
* nEq -> numero de linhas e colunas                                       *
* ----------------------------------------------------------------------- *
* Parametros de saida:                                                    *
* ----------------------------------------------------------------------- *
* n -> precondicionador                                                   *
* ----------------------------------------------------------------------- *
* OBS:                                                                    *
* ----------------------------------------------------------------------- *
***************************************************************************/
void ilu(double *RESTRICT m, double *RESTRICT a, int const nEq)
{

  int i, j, k;

/*...*/
  for(i = 0; i < nEq; i++)
    for(j = 0; j < nEq; j++)
      MAT2D(i, j, m, nEq) = MAT2D(i, j, a, nEq);
/*........................................................................*/

/*...*/
  for(k = 0; k < nEq - 1; k++){
    for(i = k + 1; i < nEq; i++){
      if (MAT2D(i, k, m, nEq) != 0.e0){
        MAT2D(i, k, m, nEq) /= MAT2D(k, k, m, nEq);  
      }
        
      for(j = k + 1; j < nEq; j++)
        if (MAT2D(i, j, m, nEq) != 0.e0)
          MAT2D(i, j, m, nEq) -= MAT2D(i, k, m, nEq) * MAT2D(k, j, m, nEq);
        
    } 
  }
/*........................................................................*/

}
/**************************************************************************/

/**************************************************************************
* DATA DE CRIACAO  : 02/10/2020                                           *
* DATA DE MODIFICAO: 00/00/0000                                           *
* ----------------------------------------------------------------------- *
* fb_ilu: forward e backward substituicao                                 *
* ----------------------------------------------------------------------- *
* Parametros de entrada:                                                  *
* ----------------------------------------------------------------------- *
* r   - vetor de indepemdente                                             *
* m   - preconsiciondador                                                 *
* x   - nao definido                                                      *
* nEq - numero de equacoes                                                *
* ----------------------------------------------------------------------- *
* Parametros de saida:                                                    *
* ----------------------------------------------------------------------- *
* x    - solucao                                                          *
* ----------------------------------------------------------------------- *
* OBS:                                                                    *
* ----------------------------------------------------------------------- *
***************************************************************************/
void fb_ilu(double *RESTRICT r, double *RESTRICT m
          , double *RESTRICT x, int const nEq)
{

  int i, j;

/*...*/
  for(i = 0; i < nEq; i++){
    x[i] = r[i];
    for(j = 0; j < i; j++)
      x[i] -=  MAT2D(i, j, m, nEq) * x[j];
  }
/*........................................................................*/

/*...*/
  for(i = nEq - 1; i > -1; i--){
    x[i] = x[i] / MAT2D(i, i, m, nEq);
    for(j = i + 1; j < nEq; j++)
      x[i] -=  MAT2D(i, j, m, nEq) * x[j] / MAT2D(i, i, m, nEq);
  } 
/*........................................................................*/

}
/**************************************************************************/




