import Pcg  as p
from scipy import io
from scipy.sparse import coo_matrix
import time as tm

def main():

# ... leitura do arquivo da matrix
  filein = 'retangulo_mec.mtx'
  aCoo  = io.mmread(filein)
  nla    = int(io.mminfo(filein)[0])
  nca    = int(io.mminfo(filein)[1])
  if nla != nca:
    print('Numero de linhas e colunas diferentes na matriz de coeficientes')
    print('Numero de linhas = {0}\nNumero de colunas'.format(nla,nca))
    exit(0)
  neq = nla

# ... leitura do arquivo do vetor de forcas
  filein = 'retangulo_mec_b.mtx'
  b      = io.mmread(filein)
  nlb    = int(io.mminfo(filein)[0])
  b      = b.reshape((nlb,))
  if nlb != neq:
    print('Numero de linhas no vertor de forcas incompativel')
    print('Numero de linhas = {0}\n'.format(nl))
    exit(0)
  neq = nlb

  print('numero de equacoes = {0}'.format(neq))
  aDense= coo_matrix(aCoo,shape=(neq,neq)).todense()

  time1 = tm.time()
  x  = p.pcgNumpy(aDense,b,preC=0)
  time1 = tm.time() - time1

  time2 = tm.time()
  x  = p.pcgNumpy(aDense,b)
  time2 = tm.time() - time2

  print('Tempo python puro = {0}\nTempo numpy       = {1}'.format(time1,time2))
  
if __name__ == '__main__':

  main()