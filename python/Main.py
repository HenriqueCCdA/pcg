from Pcg import pcg
from Pcg_numpy import pcgNumpy
from scipy import io
from scipy.sparse import coo_matrix
import time as tm
import plotLog as pLog

def main():

# ... leitura do arquivo da matrix
# filein = '../data/sist3.mtx'
  filein = '../data/retangulo_dif_4.mtx'

  aCoo  = io.mmread(filein)
  nla    = int(io.mminfo(filein)[0])
  nca    = int(io.mminfo(filein)[1])
  if nla != nca:
    print('Numero de linhas e colunas diferentes na matriz de coeficientes')
    print('Numero de linhas = {0}\nNumero de colunas'.format(nla,nca))
    exit(0)
  neq = nla

# ... leitura do arquivo do vetor de forcas
# filein = '../data/sist3_b.mtx'
  filein = '../data/retangulo_dif_4_b.mtx'

  b      = io.mmread(filein)
  nlb    = int(io.mminfo(filein)[0])
  b      = b.reshape((nlb,))
  if nlb != neq:
    print('Numero de linhas no vertor de forcas incompativel')
    print('Numero de linhas = {0}\n'.format(nlb))
    exit(0)
  neq = nlb

  print('**************************')
  print('numero de equacoes = {0}'.format(neq))
  print('**************************')
  aDense= coo_matrix(aCoo,shape=(neq,neq)).toarray()
# ... sem precondicionador
  print('Sem precondicionador')
  time1 = tm.time()
  x  = pcgNumpy(aDense,b,preC=1,nameLog='CG.txt')
#  x  = pcg(aDense,b,preC=0,nameLog='CG.txt')
  time1 = tm.time() - time1
# ......................................................................

# ... com precondicionador diagonal
  print('Precondicionador diagonal')
  time2 = tm.time()
  x  = pcgNumpy(aDense,b,preC=2,nameLog='PCG.txt')
# x  = pcg(aDense,b,preC=1,nameLog='PCG.txt')
  time2 = tm.time() - time2
# ......................................................................

  print('Tempo python puro = {0}\nTempo numpy       = {1}'.format(time1,time2))

# ... plota a convergencia
  pLog.initPlot(iCod=0)
  pLog.plot('CG.txt' )
  pLog.plot('PCG.txt')
  pLog.show()
# ......................................................................
  
if __name__ == '__main__':

  main()
