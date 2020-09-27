import numpy as np
from Precond import pre_cond_solver, pre_make, fb_dilu

def pcgNumpy(a, b, neq: int, x = None, tol: float = 1.e-14
            , maxIt: int = 10000, preC: int = 1
            ,new: bool = True, fHist: bool = True
            , nameLog = 'log.txt'):
    """
    ***************************************************************************
    * DATA DE CRIACAO  : 28/02/2019                                           *
    * DATA DE MODIFICAO: 24/09/2020                                           *
    * ----------------------------------------------------------------------- *
    * pcgNumpy : gradiente conjugado com precondicionador diagonal            *
    * ----------------------------------------------------------------------- *
    * Parametros de entrada:                                                  *
    * ----------------------------------------------------------------------- *
    * a - matriz de coefientes n x n                                          *
    * b - vetor independente                                                  *
    * neq - numero de equacoes                                               *
    * x - chute inicial                                                       *
    * tol - tolerancia do solver                                              *
    * maxIt - maximo de itracoes                                              *
    * new   - true inicialia com zeros                                        *
    * preCond - tipo do precondicionador                                      *
    *         0 - nenhum                                                      *
    *         1 - Jacobi                                                      *
    *         2 - DILU                                                        *
    * ----------------------------------------------------------------------- *
    * Parametros de saida:                                                    *
    * ----------------------------------------------------------------------- *
    * x - vetor solucao                                                       *
    * ----------------------------------------------------------------------- *
    * OBS:                                                                    *
    * ----------------------------------------------------------------------- *
    * Utilizando numpy                                                        *
    ***************************************************************************
    """
    # ... abre o arquivo de log
    if fHist :
        fileLog = open(nameLog,'w')
    # .........................................................................

    # ...
    if type(a) is type([]):
        print('Conv List a for numpy matrix')
        a = np.asarray(a)
    # ...
    if type(x) is type([]):
        print('Conv List a for numpy array')
        x = np.asarray(x)
    # ...
    if type(b) is type([]):
        print('Conv List a for numpy array')
        b = np.asarray(b)
    # .........................................................................

    # ...
    b = b.reshape((neq,1))
    # .........................................................................

    # ... precond
    m = pre_make(a, preC, neq)
    # .........................................................................

    # ...
    if new:
        x = np.zeros((neq,1),dtype=float)
    # .........................................................................

    # ...
    d0: int = float(np.dot(b.T,b))
    d: int = d0
    norm_b: int = np.sqrt(d)
    conv: int  = tol*norm_b
    # .........................................................................

    #...
    if fHist:
        fileLog.write("{0:10d},{1:.14e}\n".format(0,np.sqrt(d/d0)))
    # .........................................................................

    # ...
    # ... r0 = b - Ax0
    r = b - a@x
    # ... z0 = (M-1)r0
    z = pre_cond_solver(m, a, r, neq, preC)
    # ... p0 = r0
    p = z
    # .........................................................................

    # ... ( r(0),z(0) ) = ( r(0), (M-1)r0 )
    d = float(np.dot(r.T,z))
    # .........................................................................

    # ...
    jj = 1
    for j in range(1,maxIt):

        # ...
        #    z = np.dot(a,p)
        z = a@p
        # .....................................................................

        # ... alpha =( r(j),z(j) ) / ( Ap(j), p(j) ))
        alpha = d/float(np.dot(z.T,p))
        # .....................................................................

        # ...
        # ... x(j+1) = x(j) + alpha*p
        x += alpha*p
        # ... r(j+1) = r(j) - alpha*Ap
        r -= alpha*z
        # ... z  = (M-1)r0
        z = pre_cond_solver(m, a, r, neq, preC)
        # .....................................................................

        # ... ( r(j+1),(M-1)r(j+1) ) = ( r(j+1),z )
        di   = float(np.dot(r.T,z))
        # ... beta = ( r(j+1),(M-1)r(j+1) ) / ( r(j),r(j) )
        beta = di/d
        # .....................................................................

        # ... p(j+1) = (M-1)r(j+1) + beta*p(j) = z + beta*p(j)
        p = z + beta*p
        # .....................................................................

        # ...
        d = di
        #...
        if fHist:
            fileLog.write("{0:10d},{1:.14e}\n".format(j,np.sqrt(abs(d)/d0)))
        # .....................................................................

        # ...
        if np.sqrt(abs(d)) < conv:
            break
        # .....................................................................

        # ...
        if jj == 500:
            jj = 0
            print(j,np.sqrt(abs(d)),conv)
        jj = jj + 1
        # .....................................................................

    # .........................................................................

    # ... Energy norm: x*Kx
    z   = a@x
    xkx = float(np.dot(x.T,z))
    # .........................................................................

    # ...
    pre_name = "CG", "JCG", "DILU(0)CG", "ILU(0)CG"
    stry = f"({pre_name[preC]}) solver:\n"
    stry+= "Solver tol           = {0:e}\n"
    stry+= "Number of equations  = {1}\n"
    stry+= "Number of iterations = {2}\n"
    stry+= "xKx                  = {3:e}\n"
    print(stry.format(tol,neq,j,xkx))
    # .........................................................................

    # ... fecha o arquivo de log
    if fHist :
        fileLog.close()
    # .........................................................................

    return x
# *****************************************************************************
