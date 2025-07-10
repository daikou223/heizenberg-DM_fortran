import numpy as np
from scipy.sparse import csr_matrix, kron, eye
from functools import reduce

# スピン演算子をcsr_sparseに変換
s_x = csr_matrix(np.array([[0, 0.5], [0.5, 0]], dtype=np.complex64))
s_y = csr_matrix(np.array([[0, -0.5j], [0.5j, 0]], dtype=np.complex64))
s_z = csr_matrix(np.array([[0.5, 0], [0, -0.5]], dtype=np.complex64))
I = eye(2, format='csr', dtype=np.complex64)

siteNum = 16

def kron_list(mat_list):
    return reduce(lambda a, b: kron(a, b, format='csr'), mat_list)

# Hamiltonianの初期化（疎行列）
dim = 2 ** siteNum
hamil = csr_matrix((dim, dim), dtype=np.complex64)

# sitesとbondの定義（例として単純な1次元鎖）
sites = [i for i in range(siteNum)]
bonds = [(i, i + 1) for i in range(siteNum - 1)]

J = 1
D = 1

for mainSite, subSite in bonds:
    for op_x, op_y in [(s_x, s_x), (s_y, s_y), (s_z, s_z)]:
        mats = [I] * siteNum
        mats[mainSite] = op_x
        mats[subSite] = op_y
        hamil -= J * kron_list(mats)

    mats = [I] * siteNum
    mats[mainSite] = s_x
    mats[subSite] = s_y
    hamil += D * kron_list(mats)

    mats = [I] * siteNum
    mats[mainSite] = s_y
    mats[subSite] = s_x
    hamil -= D * kron_list(mats)

print("Hamiltonian constructed.")

# ここで直接固有値を求めるのは大変なので、疎行列用の固有値ソルバーを使うべき
from scipy.sparse.linalg import eigsh

# 固有値計算（最小の6個を求める例）
eigvals, eigvecs = eigsh(hamil, k=6, which='SA')
print("Lowest 6 eigenvalues:")
print(eigvals)
