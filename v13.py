import numpy as np
from itertools import combinations,permutations
import math
import datetime
from scipy.sparse.linalg import LinearOperator, eigsh
import os
import datetime
import time
from numba import njit
from numba.typed import List
#定数***********************************
bond_mod_J = 1
bond_mod_D = 0.1
accuracy = 1e-5 #精度
vec = None
UP_NUM = 0
inputFile = None
N = 3
maxSearch = 1
saveTimeing = 60 #second
chashFile = None
startTime = None
vectors = []

#Util関数*************************
def numberInit(Sitenum,vec):
    States = []
    idToInd = np.zeros(2**Sitenum,dtype = np.int64)
    enegys = []
    Bonds = List()
    inputFile= f"{vec[0]}_{vec[1]}_N={Sitenum}"
    flag = 1
    return States,idToInd,enegys,Bonds,inputFile,flag

def setter(SITE_NUM,inputFile,bond_mod_D):
    Bonds = List()
    stateFile = open(f'./input/N={SITE_NUM}/{inputFile}.txt',"r")
    while 1:
        line = stateFile.readline().replace("\n", "")
        if(line == ""):
            break
        [Before_bond,After_bond,mod] = list(map(int,line.split(" ")))
        Bonds.append((float(Before_bond), float(After_bond), float(mod*bond_mod_D)))
    stateFile.close()
    return Bonds

@njit
def Hamil(vecData,bond_mod_J,Bonds,UP_NUM,bond_mod_D,allComb,idToInd):
        ret_hamil = np.full(len(vecData), 0, dtype=np.complex128)
        ret_hamil += (1/4 * bond_mod_J * len(Bonds)) * vecData
        for bondInd in range(len(Bonds)):
            # 状態生成ループ
            for comb  in allComb[bondInd]:
                state = 1 << Bonds[bondInd][0]
                fripState = 1 << Bonds[bondInd][1]
                for upSite in comb:
                    state |= 1 << upSite
                    fripState |= 1 << upSite
                stateIndex = idToInd[state]
                fripStateIndex = idToInd[fripState]
                ret_hamil[stateIndex]     += -0.5 * bond_mod_J * vecData[stateIndex]
                ret_hamil[fripStateIndex] += -0.5 * bond_mod_J * vecData[fripStateIndex]
                ret_hamil[fripStateIndex] += (0.5 * bond_mod_J - 0.5j * Bonds[bondInd][2]) * vecData[stateIndex]
                ret_hamil[stateIndex]     += (0.5 * bond_mod_J + 0.5j * Bonds[bondInd][2]) * vecData[fripStateIndex]
        # if(1):
        #     path = f"./saveData/{N}_{UP_NUM}_{bond_mod_J}_{bond_mod_D}.npy"
        #     np.save(path,ret_hamil)
        #     print(f"save Site={N} in {UP_NUM} {datetime.datetime.now()}")
        return ret_hamil

def checkChash(SITE_NUM,UP_NUM,bond_mod_J,bond_mod_D):
    try:
        path = f"./saveData/{SITE_NUM}_{UP_NUM}_{bond_mod_J}_{bond_mod_D}.npy"
        v0 = np.load(path)
        return v0
    except:
        print("キャッシュファイルがありません")
        return np.array([])

def combGeneletor(N,Bonds,UP_NUM):
    ret_comb = List()
    for bond in Bonds:
        Range = set(range(N))
        Range.remove(bond[0])
        Range.remove(bond[1])
        comb_list = np.array(list(combinations(Range, UP_NUM - 1)), dtype=np.int64)
        typed_comb = List()
        for c in comb_list:
            typed_c = List()
            for x in c:
                typed_c.append(x)
        ret_comb.append(typed_comb)
    return ret_comb

#メイン処理**************************
def main():
    os.makedirs(f'./output/N={N}_aut', exist_ok=True)
    ListFile = open(f"./input/N={N}/List.txt","r")
    while 1:
        line = ListFile.readline().replace("\n", "")
        if(line == ""):
            break
        line = list(map(int,line.split(",")))
        vectors.append([(line[0],line[1]),(line[2],line[3])])
    alreadyEnegy = []
    for i in range(min(len(vectors),maxSearch)):
        print(vectors[i])
        States,idToInd,enegys,Bonds,inputFilePath,flag = numberInit(N,vectors[i])
        Bonds = setter(N,inputFilePath,bond_mod_D)
        programStart = datetime.datetime.now()
        for UP_NUM in range(1,N,1):
            DOWN_NUM = N-UP_NUM
            #逆のときは対称性より計算しなくてよい
            if(UP_NUM<=DOWN_NUM):
                startTime = datetime.datetime.now()
                States,idToInd,enegys,Bonds,inputFilePath,flag = numberInit(N,vectors[i])
                print(f'UP_NUM is {UP_NUM}')
                ALL_STATE_NUM = math.comb(N,UP_NUM)
                print(f"this hamiltonian is {ALL_STATE_NUM} squred")
                StateInd = 0
                #状態と逆変換を列挙
                for ups in combinations(range(N),UP_NUM):
                    spins = 0
                    for i in ups:
                        spins |= 1 << i
                    States.append(spins)
                    idToInd[spins] = StateInd
                    StateInd += 1
                allComb = combGeneletor(N,Bonds,UP_NUM)
                hamilton_matrix = LinearOperator(
                    shape=(ALL_STATE_NUM, ALL_STATE_NUM),
                    matvec=lambda v: Hamil(v,bond_mod_J,Bonds,UP_NUM,bond_mod_D,allComb,idToInd),
                    dtype=np.complex128
                )
                v_0 = checkChash(N,UP_NUM,bond_mod_J,bond_mod_D)
                if(len(v_0) == 0):
                    vals, vecs = eigsh(hamilton_matrix, k=1, which='SA', tol=accuracy) 
                else:
                    vals, vecs = eigsh(hamilton_matrix, v0 = v_0,k=1, which='SA', tol=accuracy) 
                print(f'MIN_ENEGY is {vals[0]}')
                enegys.append(vals[0])
                print((datetime.datetime.now()-startTime).seconds)
                print(datetime.datetime.now())
                time.sleep((datetime.datetime.now()-startTime).seconds/2)
                if(UP_NUM == 4):
                    flag = 1
                    for i in alreadyEnegy:
                        if(abs(vals[0] - i) < 10**-4):
                            flag = 0
                    if(flag == 1):
                        alreadyEnegy.append(vals[0])
                    else:
                        break
                if(flag == 1):
                    #出力ファイルの準備
                    logFile = open(f"./output/N={N}_aut/{os.path.basename(__file__)}_{inputFile}_J={bond_mod_J}_D={bond_mod_D}.txt","w")
                    logFile.write("S_Z^2,E/j\n")
                    for UP_NUM in range(1,N,1):
                        DOWN_NUM = N-UP_NUM
                        #逆のときは対称性より計算しなくてよい
                        if(UP_NUM<=DOWN_NUM and UP_NUM <= len(enegys)):
                            logFile.write(f"{(UP_NUM*1/2-DOWN_NUM*1/2)**2},{enegys[UP_NUM-1]}\n")
                    logFile.write(str(datetime.datetime.now()-programStart))
                    logFile.close()
                else:
                    print("同型のため省略")

main()