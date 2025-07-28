import numpy as np
from itertools import combinations
import math
import datetime
from scipy.sparse.linalg import LinearOperator, eigsh


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
SITE_NUM = 16
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
#ハミルトニアン操作
def Hamil(vecData):
    ret_hamil = np.zeros(len(vecData),dtype=np.complex128) 
    for ind,vecElm in enumerate(vecData):
        for bond in Bonds:
            spin1 = get_bit(States[ind],bond.beforeSite)
            spin2 = get_bit(States[ind],bond.afterSite)
            fripState = States[ind]^(1<<(SITE_NUM-1-bond.beforeSite) | 1<<(SITE_NUM-1-bond.afterSite))
            if(spin1 == spin2):
                ret_hamil[ind] += 1/4*bond_mod_J*vecElm
            elif(spin1 == 1 and spin2 == 0):
                ret_hamil[idToInd[fripState]] += (1/2*bond_mod_J-1/2j*bond.mod)*vecElm
                ret_hamil[ind] += -1/4*bond_mod_J*vecElm
            elif(spin1 == 0 and spin2 == 1):
                ret_hamil[idToInd[fripState]] += (1/2*bond_mod_J+1/2j*bond.mod)*vecElm
                ret_hamil[ind] += -1/4*bond_mod_J*vecElm
            else:
                print("この状態はないよ")
    return ret_hamil

def main():
    global SITE_NUM,States,idToInd
    #inputfileからの読み込み
    setter()
    #出力ファイルの準備
    logFile = open(f"./output/{inputFile}_J={bond_mod_J}_D={bond_mod_D}.txt","w")
    logFile.write("S_Z^2,E/j\n")
    #アップスピンとダウンスピンの数の差分DIFF_UP_DOWN_SPINに対して最低Eを取る(E(-1) = E(1))より、半分は割愛
    #DIFF_UP_DUWN_SPIN=SITE_NUMのときはすべてのスピンがそろっているので、基底になりえないから割愛
    for DIFF_UP_DOWN_SPIN in range(-SITE_NUM+2,1,2):
        startTime = datetime.datetime.now()
        #状態と逆変換を初期化
        States = []
        idToInd = {}
        print(f'DIFF_UP_DOWN_SPIN is {DIFF_UP_DOWN_SPIN}')
        ALL_STATE_NUM = math.comb(SITE_NUM,(DIFF_UP_DOWN_SPIN+SITE_NUM)//2)
        print(f"this hamiltonian is {ALL_STATE_NUM} squred")
        StateInd = 0
        #状態と逆変換を列挙
        for ups in combinations(range(SITE_NUM),(DIFF_UP_DOWN_SPIN+SITE_NUM)//2):
            spins = 0
            for i in ups:
                spins |= (1<<(SITE_NUM-1-i))
            States.append(spins)
            idToInd[spins] = StateInd
            StateInd += 1
        hamilton_matrix = LinearOperator(
            shape=(ALL_STATE_NUM, ALL_STATE_NUM),
            matvec=Hamil,
            dtype=np.complex128
        )
        vals, vecs = eigsh(hamilton_matrix, k=1, which='SA') 
        print(f'MIN_ENEGY is {vals[0]}')
        print(datetime.datetime.now()-startTime)
    logFile.close()

main()