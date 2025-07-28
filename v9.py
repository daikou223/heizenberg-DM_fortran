import numpy as np
from itertools import combinations,permutations
import math
import datetime
from scipy.sparse.linalg import LinearOperator, eigsh
import os


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
#定数***********************************
bond_mod_J = 1
bond_mod_D = 0.1
accuracy = 1e-5 #精度
Bonds = []
SITE_NUM = 21
UP_NUM = 0
inputFile = f"N={SITE_NUM}"
#変数********************************
#Stateの配列
States = []
idToInd = {}
hamil = []
rows = []
cols = []
datas = []
#関数定義******************************
def setter():
    stateFile = open("./newInput/" + inputFile+".txt","r")
    while 1:
        line = stateFile.readline().replace("\n", "")
        if(line == ""):
            break
        [Before_bond,After_bond,mod] = list(map(int,line.split(" ")))
        Bonds.append(Bond(Before_bond,After_bond,mod*bond_mod_D))
    stateFile.close()
def get_bit(num,i):
    return num >> (SITE_NUM-1-i) & 1
def disp_bit(num,length):
    bitString = bin(num)[2:]
    return "0"*(length-len(bitString))+bitString
#ハミルトニアン操作
def Hamil(vecData):
    global SITE_NUM,UP_NUM
    ret_hamil = np.full(len(vecData),0,dtype = np.complex128)
    #1,0のとき
    for bond in Bonds:
        state = 0
        fripState = 0
        Range = set(range(0,SITE_NUM))
        Range.remove(bond.beforeSite)
        Range.remove(bond.afterSite)
        ret_hamil += 1/4*bond_mod_J*vecData
        for upSites in combinations(Range,UP_NUM-1):
            state = 1 << bond.beforeSite
            fripState = 1 << bond.afterSite
            for upSite in upSites:
                state |= 1 << upSite
                fripState |= 1 << upSite
            stateIndex = idToInd[state]
            fripStateIndex = idToInd[fripState]
            ret_hamil[stateIndex] += -1/2*bond_mod_J*vecData[stateIndex]
            ret_hamil[fripStateIndex] += -1/2*bond_mod_J*vecData[fripStateIndex]
            ret_hamil[fripStateIndex] += (1/2*bond_mod_J-1/2j*bond.mod)*vecData[stateIndex]
            ret_hamil[stateIndex] += (1/2*bond_mod_J+1/2j*bond.mod)*vecData[fripStateIndex]
    return ret_hamil

def main():
    global SITE_NUM,States,idToInd,loopcounter,UP_NUM
    #inputfileからの読み込み
    setter()
    #出力ファイルの準備
    logFile = open(f"./output/{os.path.basename(__file__)}_{inputFile}_J={bond_mod_J}_D={bond_mod_D}.txt","w")
    logFile.write("S_Z^2,E/j\n")
    #アップスピンとダウンスピンの数の差分UP_NUMに対して最低Eを取る(E(-1) = E(1))より、半分は割愛
    #DIFF_UP_DUWN_SPIN=SITE_NUMのときはすべてのスピンがそろっているので、基底になりえないから割愛
    programStart = datetime.datetime.now()
    for UP_NUM in range(1,SITE_NUM,1):
        DOWN_NUM = SITE_NUM-UP_NUM
        #逆のときは対称性より計算しなくてよい
        if(UP_NUM<=DOWN_NUM):
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
            hamilton_matrix = LinearOperator(
                shape=(ALL_STATE_NUM, ALL_STATE_NUM),
                matvec=Hamil,
                dtype=np.complex128
            )
            vals, vecs = eigsh(hamilton_matrix, k=1, which='SA', tol=1e-4) 
            print(f'MIN_ENEGY is {vals[0]}')
            print(datetime.datetime.now()-startTime)
            logFile.write(f"{(UP_NUM*1/2-DOWN_NUM*1/2)**2},{vals[0]}\n")
    logFile.write(str(datetime.datetime.now()-programStart))
    logFile.close()

main()