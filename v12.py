import numpy as np
from itertools import combinations,permutations
import math
import datetime
from scipy.sparse.linalg import LinearOperator, eigsh
import os
import datetime
import time

def warraper(N,vector):
    #クラス定義***************************
    class Site():
        def __init__(self,Id,bond,mod):
            #Idはこのサイトの番号、bondはこれと結合しているサイトIdの配列(self.id<bond.id),modは定数(+上向き)
            self.id = Id
            self.bond = bond
            self.mod = mod
        def __repr__(self):
            return f"Site(id={self.id},bond = {self.bond},mod = {self.mod})"
    class Bond():
        def __init__(self,beforeSite,afterSite,bond_mod):
            #Idはこのサイトの番号、bondはこれと結合しているサイトIdの配列(self.id<bond.id),modは定数(+上向き)
            self.beforeSite = beforeSite
            self.afterSite = afterSite
            self.mod = bond_mod
        def __repr__(self):
            return f"bond(beforeSite = {self.beforeSite},afterSite = {self.afterSite},mod = {self.mod})"
    #関数定義******************************
    def setter(SITE_NUM,inputFile,bond_mod_D):
        global Bonds
        Bonds = []
        stateFile = open(f'./input/N={SITE_NUM}/{inputFile}.txt',"r")
        while 1:
            line = stateFile.readline().replace("\n", "")
            if(line == ""):
                break
            [Before_bond,After_bond,mod] = list(map(int,line.split(" ")))
            Bonds.append(Bond(Before_bond,After_bond,mod*bond_mod_D))
        stateFile.close()
    #ハミルトニアン操作
    def Hamil(vecData):
        global saveTimeing
        ret_hamil = np.full(len(vecData), 0, dtype=np.complex128)
        # 加算用バッファ
        updates_idx = []
        updates_val = []
        ret_hamil += (1/4 * bond_mod_J * len(Bonds)) * vecData
        for bond in Bonds:
            Range = set(range(SITE_NUM))
            Range.remove(bond.beforeSite)
            Range.remove(bond.afterSite)
            # 状態生成ループ
            for upSites in combinations(Range, UP_NUM - 1):
                state = 1 << bond.beforeSite
                fripState = 1 << bond.afterSite
                for upSite in upSites:
                    state |= 1 << upSite
                    fripState |= 1 << upSite
                stateIndex = idToInd[state]
                fripStateIndex = idToInd[fripState]
                # 更新をバッファに貯める
                updates_idx.append(stateIndex)
                updates_val.append(-0.5 * bond_mod_J * vecData[stateIndex])
                updates_idx.append(fripStateIndex)
                updates_val.append(-0.5 * bond_mod_J * vecData[fripStateIndex])
                updates_idx.append(fripStateIndex)
                updates_val.append((0.5 * bond_mod_J - 0.5j * bond.mod) * vecData[stateIndex])
                updates_idx.append(stateIndex)
                updates_val.append((0.5 * bond_mod_J + 0.5j * bond.mod) * vecData[fripStateIndex])
        # 一括加算（Cレベル処理で高速化）
        np.add.at(ret_hamil, updates_idx, updates_val)
        if((datetime.datetime.now()-startTime).seconds > saveTimeing):
            path = f"./saveData/{SITE_NUM}_{UP_NUM}_{bond_mod_J}_{bond_mod_D}.npy"
            np.save(path,ret_hamil)
            saveTimeing += 60
            print(f"save Site={N} in {UP_NUM} {datetime.datetime.now()}")
        return ret_hamil

    def main():
        global SITE_NUM,vec,UP_NUM,inputFile,idToInd,States,ALL_STATE_NUM,alreadyEnegy,Bonds,startTime,saveTimeing
        States = []
        idToInd = {}
        enegys = []
        SITE_NUM = N
        vec = vector
        Bonds = []
        inputFile= f"{vec[0]}_{vec[1]}_N={SITE_NUM}"
        flag = 1
        #inputfileからの読み込み
        setter(SITE_NUM,inputFile,bond_mod_D)
        #アップスピンとダウンスピンの数の差分UP_NUMに対して最低Eを取る(E(-1) = E(1))より、半分は割愛
        #DIFF_UP_DUWN_SPIN=SITE_NUMのときはすべてのスピンがそろっているので、基底になりえないから割愛
        programStart = datetime.datetime.now()
        enegys = []
        for UP_NUM in range(1,SITE_NUM,1):
            DOWN_NUM = SITE_NUM-UP_NUM
            #逆のときは対称性より計算しなくてよい
            if(UP_NUM<=DOWN_NUM):
                saveTimeing = 60
                startTime = datetime.datetime.now()
                #状態と逆変換を初期化
                States = []
                idToInd = {}
                print(f'UP_NUM is {UP_NUM}')
                ALL_STATE_NUM = math.comb(SITE_NUM,UP_NUM)
                print(f"this hamiltonian is {ALL_STATE_NUM} squred")
                StateInd = 0
                #状態と逆変換を列挙
                for ups in combinations(range(SITE_NUM),UP_NUM):
                    spins = 0
                    for i in ups:
                        spins |= 1 << i
                    States.append(spins)
                    idToInd[spins] = StateInd
                    StateInd += 1
                val = None
                val = culcMinEnegy(SITE_NUM,UP_NUM)
                # print(Hamil(np.array([1,0,0], dtype=np.complex128)))
                # print(Hamil(np.array([0,1,0], dtype=np.complex128)))
                # print(Hamil(np.array([0,0,1], dtype=np.complex128)))
                if(val is not None):
                    enegys.append(val)
                    print((datetime.datetime.now()-startTime).seconds)
                    print(datetime.datetime.now())
                    time.sleep((datetime.datetime.now()-startTime).seconds/2)
                    if(UP_NUM == 4):
                        flag = 1
                        for i in alreadyEnegy:
                            if(abs(val - i) < 10**-4):
                                flag = 0
                        if(flag == 1):
                            alreadyEnegy.append(val)
                        else:
                            break
        if(flag == 1):
            #出力ファイルの準備
            logFile = open(f"./output/N={SITE_NUM}_aut/{os.path.basename(__file__)}_{inputFile}_J={bond_mod_J}_D={bond_mod_D}.txt","w")
            logFile.write("S_Z^2,E/j\n")
            for UP_NUM in range(1,SITE_NUM,1):
                DOWN_NUM = SITE_NUM-UP_NUM
                #逆のときは対称性より計算しなくてよい
                if(UP_NUM<=DOWN_NUM and UP_NUM <= len(enegys)):
                    logFile.write(f"{(UP_NUM*1/2-DOWN_NUM*1/2)**2},{enegys[UP_NUM-1]}\n")
            logFile.write(str(datetime.datetime.now()-programStart))
            logFile.close()
        else:
            print("同型のため省略")

    def culcMinEnegy(SITE_NUM,UP_NUM):
        global saveTimeing
        hamilton_matrix = LinearOperator(
            shape=(ALL_STATE_NUM, ALL_STATE_NUM),
            matvec=Hamil,
            dtype=np.complex128
        )
        v_0 = checkChash(SITE_NUM,UP_NUM)
        if(len(v_0) == 0):
            vals, vecs = eigsh(hamilton_matrix, k=1, which='SA', tol=accuracy) 
        else:
            vals, vecs = eigsh(hamilton_matrix, v0 = v_0,k=1, which='SA', tol=accuracy) 
        print(saveTimeing)
        if((datetime.datetime.now()-startTime).seconds > saveTimeing):
            path = f'./saveData/{SITE_NUM}_{UP_NUM}_{bond_mod_J}_{bond_mod_D}.npy'
            np.save(path,vecs[0])
            saveTimeing += 60
            print(f"save Site={N} in {UP_NUM} {datetime.datetime.now()}")
        print(f'MIN_ENEGY is {vals[0]}')
        return vals[0]
    
    def checkChash(SITE_NUM,UP_NUM):
        try:
            path = f"./saveData/{SITE_NUM}_{UP_NUM}_{bond_mod_J}_{bond_mod_D}.npy"
            v0 = np.load(path)
            return v0
        except:
            print("キャッシュファイルがありません")
            return np.array([])
    main()


#定数***********************************
bond_mod_J = 1
bond_mod_D = 0.1
accuracy = 1e-5 #精度
SITE_NUM = None
vec = None
UP_NUM = 0
inputFile = None
N = 3
maxSearch = 1
saveTimeing = 60 #second
chashFile = None
startTime = None
#変数********************************
#Stateの配列
States = []
Bonds = []
idToInd = {}
hamil = []
rows = []
cols = []
datas = []
enegys = []
before_enegy = None
vectors = []
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
    warraper(N,vectors[i])
