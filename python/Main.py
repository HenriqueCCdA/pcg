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
#  file_in_a = '../data/cilindro1_1_1.mtx'
#  file_in_b = '../data/cilindro1_1_1_b.mtx'
#  file_in_a = '../data/sist3.mtx'
#  file_in_b = '../data/sist3_b.mtx'

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

  # convertendo para uma matriz cheia
  aDense= coo_matrix(aCoo,shape=(neq,neq)).toarray()

  return aDense, b, neq

def main():

  a, b, neq = read_mm()

  # ...
  time = []
  list_name = ('CG', 'JCG', 'DILUCG', 'ILUCG')

  for pre, name in zip([0, 1, 2, 3], list_name):
    print(name)
    time1 = tm.time()
    x  = pcgNumpy(a, b, neq, preC = pre, tol = 1e-11,  nameLog =  name +'.txt')
    time.append(tm.time() - time1)
# ......................................................................

  for name, t in zip(list_name, time):
    print(f'{name:8} = {t:.6f}')

# ... plota a convergencia
  pLog.initPlot(iCod = 0)
  for name in list_name:
    pLog.plot(name+'.txt' )
  pLog.show()
# ......................................................................
  
if __name__ == '__main__':

  main()
