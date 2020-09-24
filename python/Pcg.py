import Blas as bl
import math as ma

# *****************************************************************************
def preDiag(a,neq):
    """
    *********************************************************************
    * DATA DE CRIACAO  : 28/02/2019                                     *
    * DATA DE MODIFICAO: 00/00/0000                                     *
    * ----------------------------------------------------------------- *
    * preDiag: monta o precondicionador diagonal (Jacobi)               *
    * ----------------------------------------------------------------- *
    * Parametros de entrada:                                            *
    * ----------------------------------------------------------------- *
    * a - matriz de coefientes n x n                                    *
    * neq - numero de linhas e colunas                                  *
    * ----------------------------------------------------------------- *
    * Parametros de saida:                                              *
    * ----------------------------------------------------------------- *
    * m    - retorna o inverso do coefiente da diagonal                 *
    * ----------------------------------------------------------------- *
    * OBS:                                                              *
    * ----------------------------------------------------------------- *
    * m e um vetor e a uma matriz quadrade cheia                        *
    *********************************************************************
      """
    m = [0.e0]*neq

    for i in range(0,neq):
        m[i] = 1.0/a[i][i]

    return m
# *****************************************************************************
#
# *****************************************************************************
def pcg(a,b,x=None,tol=1.e-14,maxIt=10000, preC=1
       ,new=True,fHist=True,nameLog='log.Txt'):
    """
    *********************************************************************
    * DATA DE CRIACAO  : 28/02/2019                                     *
    * DATA DE MODIFICAO: 24/09/2020                                     *
    * ----------------------------------------------------------------- *
    * pcg : gradiente conjugado com precondicionador diagonal           *
    * ----------------------------------------------------------------- *
    * Parametros de entrada:                                            *
    * ----------------------------------------------------------------- *
    * a - matriz de coefientes n x n                                    *
    * b - vetor independente                                            *
    * x - chute inicial                                                 *
    * tol - tolerancia do solver                                        *
    * maxIt - maximo de itracoes                                        *
    * new   - true inicialia com zeros                                  *
    * ----------------------------------------------------------------- *
    * Parametros de saida:                                              *
    * ----------------------------------------------------------------- *
    * x - vetor solucao                                                 *
    * ----------------------------------------------------------------- *
    * OBS:                                                              *
    * ----------------------------------------------------------------- *
    * Utilizando python puro                                            *
    *********************************************************************
    """
# ... abre o arquivo de log
    if fHist :
        fileLog = open(nameLog,'w')
# .............................................................................

# ...
    neq = len(b)
# .............................................................................

# ...
    a = a.tolist()
    b = b.tolist()
# .............................................................................

# ... precond Jacobi
    if preC == 1:
        m = preDiag(a,neq)
    elif preC == 0:
        m = [1.0]*neq
# .............................................................................

# ...
    z = [0.0]*neq
    r = [0.0]*neq
    p = [0.0]*neq
# .............................................................................

# ...
    if x == None:
        x = neq*[0.e0]
    else:
        if new:
            for i in range(0,neq):
                x[i] = 0.e0
# .............................................................................

# ...
    for i in range(0,neq):
        z[i] = b[i]*m[i]
# .............................................................................

# ...
    d0 = d = bl.dot1(b,b)
    norm_b = ma.sqrt(d)
    conv   = tol*norm_b
# .............................................................................

# ... Ax0
    z = bl.matvec1(a,x)
# .............................................................................

#...
    if fHist:
        fileLog.write("{0:10d},{1:.14e}\n".format(0,ma.sqrt(d/d0)))
# .............................................................................

# ...
    for i in range(0,neq):
# ... r0 = b - Ax0
        r[i] = b[i] - z[i]
# ... z0 = (M-1)r0
        z[i] = r[i]*m[i]
# ... p0 = r0
        p[i] = z[i]
# .............................................................................

# ... ( r(0),z(0) ) = ( r(0), (M-1)r0 )
    d = bl.dot1(r,z)
# .............................................................................

# ...
    jj = 1
    for j in range(1,maxIt):
# ... z = Ap
        z = bl.matvec1(a,p)
# .............................................................................

# ... alpha =( r(j),z(j) ) / ( Ap(j), p(j) ))
        alpha = d/bl.dot1(z,p)
# .............................................................................

# ...
        for i in range(0,neq):
# ... x(j+1) = x(j) + alpha*p
            x[i] = x[i] + alpha*p[i]
# ... r(j+1) = r(j) - alpha*Ap
            r[i] = r[i] - alpha*z[i]
# ... z  = (M-1)r0
            z[i] = r[i]*m[i]
# .............................................................................

# ... ( r(j+1),(M-1)r(j+1) ) = ( r(j+1),z )
        di   = bl.dot1(r,z)
# ... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) ) 
        beta = di/d
# .............................................................................

# ...
        for i in range(0,neq):
# ... p(j+1) = (M-1)r(j+1) + beta*p(j) = z + beta*p(j)
            p[i] = z[i] + beta*p[i]
# .............................................................................

# ...    
        d = di
#...
        if fHist:
            fileLog.write("{0:10d},{1:.14e}\n".format(j,ma.sqrt(d/d0)))
# .............................................................................

# ...
        if ma.sqrt(d) < conv:
            break
# .............................................................................

# ...
        if jj == 500 :
            jj = 0
            print(j,ma.sqrt(d),conv)
        jj = jj + 1
# .............................................................................

# .............................................................................

# ... Energy norm: x*Kx
    z = bl.matvec2(a,x)
    xkx = bl.dot1(x,z)
# .............................................................................

# ...
    stry = "(PCG) solver:\n"
    stry+= "Solver tol           = {0:e}\n"
    stry+= "Number of equations  = {1}\n"
    stry+= "Number of iterations = {2}\n"
    stry+= "xKx                  = {3:e}\n"
    print(stry.format(tol,neq,j,xkx))
# .............................................................................
        
    return x
# *****************************************************************************
