from Pcg import pcg
from Pcg_numpy import pcgNumpy
from scipy import io
from scipy.sparse import coo_matrix
import time as tm
import plotLog as pLog

def read_mm():

  # ... leitura do arquivo da matrix
  file_in_a = '../data/retangulo_dif_4.mtx'
  file_in_b = '../data/retangulo_dif_4_b.mtx'

  aCoo  = io.mmread(file_in_a)
  nla    = int(io.mminfo(file_in_a)[0])
  nca    = int(io.mminfo(file_in_a)[1])
  if nla != nca:
    print('Numero de linhas e colunas diferentes na matriz de coeficientes')
    print('Numero de linhas = {0}\nNumero de colunas'.format(nla,nca))
    exit(0)
  neq = nla

  # ... leitura do arquivo do vetor de forcas
  b      = io.mmread(file_in_b)
  nlb    = int(io.mminfo(file_in_b)[0])
  b      = b.reshape((nlb,))
  if nlb != neq:
    print('Numero de linhas no vertor de forcas incompativel')
    print('Numero de linhas = {0}\n'.format(nlb))
    exit(0)

  print('**************************')
  print('numero de equacoes = {0}'.format(neq))
  print('**************************')
  aDense= coo_matrix(aCoo,shape=(neq,neq)).toarray()

  return aDense, b, neq

def main():

  a, b, neq = read_mm()

# ... sem precondicionador
  print('Sem precondicionador')
  time1 = tm.time()
  x  = pcgNumpy(a, b, neq, preC = 1, nameLog = 'CG.txt')
#  x  = pcg(aDense,b,preC=0,nameLog='CG.txt')
  time1 = tm.time() - time1
# ......................................................................

# ... com precondicionador diagonal
  print('Precondicionador diagonal')
  time2 = tm.time()
  x  = pcgNumpy(a, b, neq, preC = 2, nameLog = 'PCG.txt')
# x  = pcg(aDense,b,preC=1,nameLog='PCG.txt')
  time2 = tm.time() - time2
# ......................................................................

  print('Tempo python puro = {0}\nTempo numpy       = {1}'.format(time1,time2))

# ... plota a convergencia
  pLog.initPlot(iCod = 0)
  pLog.plot('CG.txt' )
  pLog.plot('PCG.txt')
  pLog.show()
# ......................................................................
  
if __name__ == '__main__':

  main()
