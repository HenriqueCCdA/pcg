#include"../include/Solv.h"
#include"../include/HccaTime.h"

/**********************************************************************
 * PCG  : metodo do gradiente conjugado com precondiconador diagonal  *
 * (M-1Ax=M-1b) (matriz simentrica)                                   *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 *  neq -> numero de equacoes                                         *
 *  nad -> numero de elementos nao nulos fora da diagonal             *
 *  ia  -> estrutura de dados para matriz esparsa A                   *
 *  ja  -> estrutura de dados para matriz esparsa A                   *
 *  al  -> parte inferior da matriz A                                 *
 *  ad  -> diagnal da matriz A                                        *
 *  au  -> parte superior da matriz A                                 *
 *   p  -> precondiconador diagonal                                   *
 *   b  -> vetor b (Ax=b)                                             *
 *   x  -> vetor de solucao                                           *
 *   z  -> vetor auxiliar                                             *
 *   r  -> vetor auxiliar                                             *
 * newX -> vetor inicial iniciado com zero                            *
 * fLog -> arquivo de log do solver                                   *
 * log  -> log de arquivo (true|false)                                *
 * tol  -> tolerancia do solver                                       *
 *maxIt -> numero maximo de iteracoes                                 *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * x[]-> atualizado                                                   *  
 * b[]-> alterado                                                     *
 * ad,al,au-> inalterado                                              *
 * -------------------------------------------------------------------*
 *********************************************************************/
void pcg(int const nEq      , double *a 
        , double *b          , double *x 
        , double const tol   , char preC
        , unsigned int maxIt , bool newX          
        , FILE* fLog         , bool log
        , void(*matvec)()    , double(*dot)())
{
  int i,j;
  int k=0;
  char transp = 1;
  double alpha, beta, d, di, conv;
  double normB, normBm, xKx, normRm, normR, normX;
  double timei,timef;
  double *z, *m, *r, *p;

  timei = getTimeC();

/*... alocando memoria*/
  z = (double *) malloc(sizeof(double)*nEq);
  ERRO_MALLOC(z, "z",__LINE__, __FILE__, __func__);
  r = (double *) malloc(sizeof(double)*nEq);
  ERRO_MALLOC(r, "r",__LINE__, __FILE__, __func__);
  p = (double *) malloc(sizeof(double)*nEq);
  ERRO_MALLOC(p, "p",__LINE__, __FILE__, __func__);
/*...................................................................*/

/*... gerando o precondicionador*/
  m = preMake( a, nEq, preC);
/*...................................................................*/

/*... chute inicial*/
	if (newX)
		for (i = 0; i < nEq; i++)
			x[i] = 0.e0;
/*...................................................................*/

/*... conv = tol * |(M-1)b|m = tol(b,M(-1)b) */
	for (i = 0; i < nEq; i++)
    z[i] = b[i] * m[i];

  d      = dot(b, z, nEq);
	normBm = sqrt(fabs(d));
	conv   = tol * normBm;
/*...................................................................*/
  
/*... ||b||*/
  normB = sqrt(dot(b, b, nEq));
/*...................................................................*/
 
/*... Ax0*/
  matvec( a, x, z, nEq, nEq, transp);
  
/*... r0 = b - Ax0*/
  for(i = 0; i < nEq; i++)   {
    r[i] = b[i] - z[i];
  }

/* ... z0 = (M-1)r0*/
  preCondSolver(m, a, r, z, nEq, preC);


/*...p0 = r0*/
 for(i = 0; i < nEq; i++)   {
    p[i] = z[i];
  }
/*...................................................................*/

/* ... ( r(0),z(0) ) = ( r(0), (M-1)r0 )*/
  d     = dot(r, z, nEq);
/*...................................................................*/

/* ... iteracao principal*/   
  for(j = 0; j < maxIt; j++)   {
    
/*... z = Ap(j)*/
    matvec(a, p, z, nEq, nEq, transp);
/*...................................................................*/

/*... alpha =( r(j),z(j) ) / ( Ap(j), p(j) )) */
    alpha = d / dot(p, z, nEq);
/*...................................................................*/
    
/*...*/
    for(i = 0; i < nEq; i++)   {
/* ... x(j+1) = x(j) + alpha*p */
      x[i] +=  alpha * p[i];
/* ... r(j+1) = r(j) - alpha*Ap*/
      r[i] -=  alpha * z[i];
    }
/*...................................................................*/

/*... z  = (M-1)r0 */
    preCondSolver(m, a, r, z, nEq, preC);
/*...................................................................*/

/* ... (r(j + 1), (M - 1)r(j + 1)) = (r(j + 1), z)*/
		di = dot(r, z, nEq);
/*... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) */
		beta = di / d;
/*...................................................................*/

/*... p(j+1) = (M-1)r(j+1) + beta*p(j) = z + beta*p(j)*/
    for(i = 0; i < nEq; i++)
      p[i] = z[i] + beta * p[i];
/*...................................................................*/   
 
/*...*/
    d = di;
    if (sqrt(fabs(d)) < conv) break;
/*....................................................................*/

/*...*/
    if( k == 500)
    { 
      printf("it: %d %e %e\n",j,sqrt(fabs(d)),conv);
      k = 0; 
    }
    k++;
/*....................................................................*/    
    
/*...*/
    if(log)
      fprintf(fLog,"%9d %e\n", j+1, sqrt(fabs(d)) / normBm);
/*....................................................................*/    
  }
/*....................................................................*/    
  
/*... Energy norm:  x*Kx*/
  matvec(a , x, z, nEq, nEq, transp);
/*... norma de energia = xTAx* */
  xKx = dot( x, z, nEq);
/*...................................................................*/

/*... norm - 2 = || x ||*/
	normX = sqrt(dot(x, x, nEq));
/*...................................................................*/
  
/*... r = M(-1)(b - Ax) (calculo do residuo explicito)*/
  for(i = 0; i < nEq; i++)
    r[i]  = b[i] - z[i];

  preCondSolver(m, a, r, z, nEq, preC);
/*...................................................................*/
  normRm = dot( r, z, nEq);
  normRm = sqrt(fabs(normRm));
  normR  = sqrt(dot( r, r, nEq));
/*...................................................................*/

  timef = getTimeC() - timei;   

/*...*/
  printf(" (CG) solver:\n"
         "\tIterarions    =      %20d\n"
         "\tEquations     =      %20d\n"
         "\tSolver conv   =      %20.6e\n"
         "\tSolver tol    =      %20.6e\n"
	       "\tx*Kx          =      %20.6e\n"
         "\t|| x ||       =      %20.6e\n"
         "\t|| b - Ax ||  =      %20.6e\n"
         "\t|| b - Ax ||m =      %20.6e\n"
         "\t||b||         =      %20.6e\n"
	       "\tCPU time(s)   =      %20.5lf\n" 
	       , j + 1, nEq,conv, tol, xKx, normX, normR, normRm, normB, timef);
  
  if(j == maxIt)
  { 
    printf(" *** WARNING: no convergende reached !\n");
    printf("MAXIT = %d \n",maxIt);
    exit(EXIT_FAILURE);
  }
/*...................................................................*/

/*... liberando memoria*/
  free(z);
  free(m);
  free(r);
  free(p);
/*....................................................................*/
}
/**********************************************************************/