#include"../include/ReadFile.h"

/**********************************************************************
 * Data de Criacao  : 28/09/2020                                      *
 * Data de Modificao: 00/00/0000                                      *
 * -------------------------------------------------------------------*
 * denseMatrix : transforma a matriz do formato COO para o cheio      *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 * a     -> nao definido                                              *
 * aCoo  -> matriz no formato COO                                     *
 * ia    -> linha do coeficiente a(k)                                 * 
 * ja    -> coluna do coeficiente a(k)                                *
 * nnz   -> numero total de numero nao nulos                          *
 * nEq   -> numero de equacoes                                        *
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * a -> matriz cheia                                                  *
 * -------------------------------------------------------------------*
 * OBS:                                                               *
 * -------------------------------------------------------------------*
 **********************************************************************/
void denseMatrix(double *a       , double const *const aCoo
               , int const * const ia, int const *const ja
               , int const nnz   , int const nEq)
{
  int i, j, k;
  
  for(k = 0; k < nEq*nEq; k++)
    a[k] = 0.0;

  for(k = 0; k < nnz; k++)
  {
    i = ia[k] - 1;
    j = ja[k] - 1;
    a[i*nEq + j] = aCoo[k];
    a[j*nEq + i] = aCoo[k];
  }
}
/**************************************************************************/



/**********************************************************************
 * Data de Criacao  : 27/09/2020                                      *
 * Data de Modificao: 00/00/0000                                      *
 * -------------------------------------------------------------------*
 * readB : le o valores do vetor B                                    *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 * b     -> vetor b (alocado)                                         *
 * fileIn-> arquivo sendo lido                                        *
 * n     -> numero de entradas lidas                                  * 
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * b -> vetor atualizado                                              *
 * -------------------------------------------------------------------*
 **********************************************************************/
int readB(double *b,FILE *fileIn,int const n){
  int i;
  int erro;
  double value;
  for(i=0;i<n;i++){
    erro = fscanf(fileIn,"%lf",&value);
    if(erro != 1) return -1;
    b[i] = value;
  }
  return 0; 
} 
/**************************************************************************/

/**********************************************************************
 * Data de Criacao  : 27/09/2020                                      *
 * Data de Modificao: 00/00/0000                                      *
 * -------------------------------------------------------------------*
 * readSystem: le o sistema de um arquivo MatrixMarket                *
 * -------------------------------------------------------------------*
 * Parametros de Entrada:                                             *
 * -------------------------------------------------------------------*
 * nEq -> nao definido                                                * 
 * -------------------------------------------------------------------*
 * Parametros de Saida:                                               *
 * -------------------------------------------------------------------*
 * nEq -> numero de equacoes                                          *
 * -------------------------------------------------------------------*
 **********************************************************************/
void readSystem(double *a[], double *b[], double *x[], int *nEq) {
  MM_typecode matCode;
  FILE *fileInA = NULL, *fileInB = NULL;
  int nLin, nCol, nnz;
  int *ia, *ja;
  double *aCoo;

  /*... open file matrix A*/
//fileInA = fopen("retangulo_mec.mtx","r");
  fileInA = fopen("cilindro1_1_1.mtx","r");
//fileInA = fopen("retangulo_dif_4.mtx","r");
  if(fileInA == NULL){ 
    printf("Erro na abertura da arquivo %s.\n","teste");
    exit(EXIT_FAILURE);
  }
 
//fileInB = fopen("retangulo_mec_b.mtx","r"); 
  fileInB = fopen("cilindro1_1_1_b.mtx","r");
//fileInB = fopen("retangulo_dif_4_b.mtx","r");
  if(fileInB == NULL){ 
    printf("Erro na abertura da arquivo %s.\n","teste");
    exit(EXIT_FAILURE);
  }

  /* checa o arquivo do matriz a*/
  mm_read_banner(fileInA,&matCode);
  mm_read_mtx_crd_size(fileInA,&nLin,&nCol,&nnz);
  
  if(! mm_is_symmetric(matCode))
  {
    printf("A matriz não é simetrica\n");
    exit(EXIT_FAILURE);
  }
/*...................................................................*/

  *nEq = nLin;


/*... checa o arquivo do vetor b*/
  mm_read_banner(fileInB,&matCode);
  mm_read_mtx_array_size(fileInB,&nLin,&nCol);

  if (*nEq != nLin)
  {
    printf("Numero de equacao da matriz A e b incompativel.\n");
    exit(EXIT_FAILURE);
  }
/*...................................................................*/


/*... alocacao da memoria*/
/*... vetor de linhas*/
  ia = (int*) malloc(sizeof(int)*nnz);
  if( ia == NULL ){
    printf("Erro na alocaco do vetor ia.\n");
    exit(EXIT_FAILURE);
  }
/*... vetor de colunas*/
  ja = (int*) malloc(sizeof(int)*nnz);
  if( ja == NULL ){
    printf("Erro na alocaco do vetor ja.\n");
    exit(EXIT_FAILURE);
  }
/*... vetor dos valores dos coeficientes*/
  aCoo  = (double*) malloc(sizeof(double)*nnz);
  if(  aCoo == NULL ){
    printf("Erro na alocaco do vetor  aCoo.\n");
    exit(EXIT_FAILURE);
  }
/*... vetor dos valores dos coeficientes*/
  *b  = (double*) malloc(sizeof(double)*(*nEq));
  if(  *b == NULL ){
    printf("Erro na alocaco do vetor  b.\n");
    exit(EXIT_FAILURE);
  }
/*... vetor dos valores dos coeficientes*/
  *a  = (double*) malloc(sizeof(double)*(*nEq)*(*nEq));
  if(  *a == NULL ){
    printf("Erro na alocaco do vetor  a.\n");
    exit(EXIT_FAILURE);
  }
/*... vetor de solucao*/
  *x  = (double*) malloc(sizeof(double)*(*nEq));
  if( *x == NULL ){
    printf("Erro na alocaco do vetor  x.\n");
    exit(EXIT_FAILURE);
  }
/*.....................................................................*/


/*... leitura da matriz*/
  mm_read_mtx_crd_data(fileInA,nLin,nCol,nnz,ia,ja,aCoo,matCode);
/*...................................................................*/

/*... leitura do vetor B*/
  readB(*b,fileInB,*nEq);
/*...................................................................*/

/*... transformando a matriz para o formato denso*/
  denseMatrix(*a, aCoo, ia, ja, nnz, *nEq);
/*...................................................................*/

/*... liberando memoria*/
  free(aCoo);
  free(ia);
  free(ja);
/*...................................................................*/

/*... fechando os arquivos*/
  fclose(fileInA);
  fclose(fileInB);
/*...................................................................*/
}
/**************************************************************************/