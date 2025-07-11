import numpy as np
from scipy.sparse import csr_matrix, csc_matrix, identity
from scipy.sparse.linalg import inv, eigs
from scipy.sparse.linalg import spsolve

def koyuuti(hamil,n):

    I = identity(n, format='csc')
    x = np.array([1+1j]*n)

    row_sums = np.array(hamil.sum(axis=1)).flatten()
    sigma = max(map(abs,row_sums))
    hamil = -hamil + (sigma+1)*I
    i = 0
    while 1:
        y = x.copy()
        x = hamil@x
        x = x / np.linalg.norm(x)
        if (np.linalg.norm(x - y) < 1e-6 or np.linalg.norm(x + y) < 1e-6) and i>30:
             break
        i += 1
    val = ((hamil @ x)[0]/x[0] - sigma - 1)*-1
    return val
