import numpy as np
from functools import reduce

np.set_printoptions(precision=3, suppress=True)
#定数**************************
J = 1
D = 1
s_x = [[0,1/2],[1/2,0]]
s_y = [[0,-1j/2],[1j/2,0]]
s_z = [[1/2,0],[0,-1/2]]
I = [[1,0],[0,1]]
#クラス定義***************************
class Site():
    def __init__(self,Id,bond):
        self.id = Id
        self.bond = bond
#系の設定**********************************
sites = [
    Site(0,[]),
    Site(1,[0])
]
siteNum = len(sites)
print("complete setting system")
#メイン処理******************************
hamil = np.zeros((2**siteNum, 2**siteNum), dtype=np.complex64)
for mainSite in range(siteNum):
    for subSite in sites[mainSite].bond:
        #内積x
        tyokuList = [I for _ in range(siteNum)]
        tyokuList[mainSite] = s_x
        tyokuList[subSite] = s_x
        hamil += -J*reduce(np.kron,tyokuList)
        print(hamil)
        #内積y
        tyokuList = [I for _ in range(siteNum)]
        tyokuList[mainSite] = s_y
        tyokuList[subSite] = s_y
        y_hamil = -J*reduce(np.kron,tyokuList)
        hamil += y_hamil
        print(y_hamil)
        #内積z
        tyokuList = [I for _ in range(siteNum)]
        tyokuList[mainSite] = s_z
        tyokuList[subSite] = s_z
        hamil += -J*reduce(np.kron,tyokuList)
        #外積
        tyokuList = [I for _ in range(siteNum)]
        tyokuList[mainSite] = s_x
        tyokuList[subSite] = s_y
        hamil += D*reduce(np.kron,tyokuList)
        tyokuList = [I for _ in range(siteNum)]
        tyokuList[mainSite] = s_y
        tyokuList[subSite] = s_x
        hamil += -D*reduce(np.kron,tyokuList)
print("hamiltonian was created")
print(hamil)
[koyuuti,vec] = np.linalg.eig(hamil)
koyuuti = sorted(koyuuti,key = lambda i:i.real)
for i in range(len(koyuuti)):
    if(koyuuti[i].imag > 1):
        print("imag part exact")
dis = 0
print(len(koyuuti)," exact")
while len(koyuuti)>dis:
    print(dis,koyuuti[dis].real)
    dis += 10
print(len(koyuuti)-1,koyuuti[-1].real)