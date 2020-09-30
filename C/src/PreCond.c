#include"../include/PreCond.h"

/***************************************************************************
* DATA DE CRIACAO  : 29/09/2019                                           *
* DATA DE MODIFICAO: 00/00/0000                                           *
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
double* preMake(double *a, int const nEq, char preC)
{
  int i;
  double *m = NULL;

/*... nenhum*/
  if(preC == 0){
    m = (double *) malloc(sizeof(double)*nEq);
    ERRO_MALLOC(m, "m",__LINE__, __FILE__, __func__);
    for(i = 0; i < nEq; i++)
      m[i] = 1.e0;
  }

/*... diagonal*/
  if(preC == 1){
    m = (double *) malloc(sizeof(double)*nEq);
    ERRO_MALLOC(m, "m",__LINE__, __FILE__, __func__);
    for(i = 0; i < nEq; i++)
      m[i] = 1.e0/a[i*nEq + i];
  }

  return m;

}
/**************************************************************************/

/**************************************************************************
* DATA DE CRIACAO  : 30/09/2019                                           *
* DATA DE MODIFICAO: 00/00/0000                                           *
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
void preCondSolver(double *m     , double *a
                 , double *r    , double *x
                 , int const nEq, char preC)
{

  int i;

/*... nenhum*/
  if(preC == 0){
    for(i = 0; i < nEq; i++)
      x[i] = r[i];
  }

/*... diagonal*/
  if(preC == 1){
    for(i = 0; i < nEq; i++)
      x[i] = r[i]*m[i];
  }


}
/**************************************************************************/