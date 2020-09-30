import numpy as np
from Precond import pre_cond_solver, pre_make, fb_dilu

def pcgNumpy(a, b, neq: int, x = None, tol: float = 1.e-14
            , maxIt: int = 10000, preC: int = 1
            ,new: bool = True, fHist: bool = True
            , nameLog = 'log.txt'):
    """
    ***************************************************************************
    * DATA DE CRIACAO  : 28/02/2019                                           *
    * DATA DE MODIFICAO: 29/09/2020                                           *
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
    *         4 - ILU                                                         *
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

    # ... conv = tol * |(M-1)b|m = tol(b,M(-1)b)/
    z = pre_cond_solver(m, a, b, neq, preC)
    d       : float = float(np.dot(b.T,z))
    norm_b_m: float = np.sqrt(d)
    conv: int  = tol*norm_b_m
    # .........................................................................

    # ... ||b||
    norm_b: int = np.sqrt(float(np.dot(b.T, b)))
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

    #...
    if fHist:
        fileLog.write("{0:10d},{1:.14e}\n".format(0,np.sqrt(d/norm_b_m)))
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
        di   = float(np.dot(r.T, z))
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
            fileLog.write("{0:10d},{1:.14e}\n".format(j,np.sqrt(abs(d)/norm_b_m)))
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
    xkx: float = float(np.dot(x.T,z))
    # .........................................................................

    # ...
    norm_x: float = np.sqrt(float(np.dot(x.T, x)))
    # .........................................................................

    # ... r =M(-1)(b - Ax)
    r = b - a@x
    z = pre_cond_solver(m, a, r, neq, preC)
    norm_r: float = np.sqrt(float(np.dot(r.T, r)))
    norm_r_m: float = np.sqrt(float(np.dot(r.T, z)))
    # .........................................................................

    # ...
    pre_name = "CG", "JCG", "DILU(0)CG", "ILU(0)CG"
    stry = f"({pre_name[preC]}) solver:\n"
    stry+= "Number of iterations = {0}\n"
    stry+= "Number of equations  = {1}\n"
    stry+= "Solver conv          = {2:e}\n"
    stry+= "Solver tol           = {3:e}\n"
    stry+= "xKx                  = {4:e}\n"
    stry+= "|| x ||              = {5:e}\n"
    stry+= "|| b - Ax ||         = {6:e}\n"
    stry+= "|| b - Ax ||m        = {7:e}\n"
    stry+= "|| b ||              = {8:e}\n"
    print(stry.format(j, neq, conv, tol, xkx, norm_x, norm_r, norm_r_m, norm_b))
    # .........................................................................

    # ... fecha o arquivo de log
    if fHist :
        fileLog.close()
    # .........................................................................

    return x
# *****************************************************************************
