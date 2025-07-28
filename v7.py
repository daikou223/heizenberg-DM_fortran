import numpy as np
from itertools import combinations
import math
import datetime
from scipy.sparse.linalg import LinearOperator, eigsh
import sys


#クラス定義***************************
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
SITE_NUM = 4
inputFile = f"N={SITE_NUM}"
#変数********************************
#Stateの配列
States = []
idToInd = {}
hamil = []
rows = []
cols = []
datas = []
chach = {}
#関数定義******************************
def setter():
    stateFile = open("./newInput/" + inputFile+".txt","r")
    while 1:
        line = stateFile.readline().replace("\n", "")
        if(line == ""):
            break
        [Before_bond,After_bond,mod] = list(map(int,line.split(" ")))
        Bonds.append(Bond(Before_bond,After_bond,mod))
    stateFile.close()
def get_bit(num,i):
    return num >> (SITE_NUM-1-i) & 1
#ハミルトニアン操作
def Hamil(vecData):
    ret_hamil = np.zeros(len(vecData),dtype=np.complex128) 
    for ind,vecElm in enumerate(vecData):
        #キャッシュには{遷移先状態1:[係数実部の4倍、係数虚部の4倍],遷移先状態2:[実部4倍、虚部4倍],...}を保存しておく
        if(ind in chach):
            for ret_ind,mod_id in chach[ind].items():
                ret_hamil[ret_ind] += (mod_id[0]*bond_mod_J/4+mod_id[1]*bond_mod_D/4*1j)*vecElm
        else:
            chachList = {}
            for bond in Bonds:
                spin1 = get_bit(States[ind],bond.beforeSite)
                spin2 = get_bit(States[ind],bond.afterSite)
                if(ind not in chachList):
                    chachList[ind] = [0,0]
                if(spin1 == spin2):
                    ret_hamil[ind] += 1/4*bond_mod_J*vecElm
                    chachList[ind][0] += 1
                else:
                    fripState = States[ind]^(1<<(SITE_NUM-1-bond.beforeSite) | 1<<(SITE_NUM-1-bond.afterSite))
                    fripInd = idToInd[fripState]
                    if(fripInd  not in chachList):
                        chachList[fripInd] = [0,0]
                    if(spin1 == 1 and spin2 == 0):
                        ret_hamil[fripInd] += (1/2*bond_mod_J-1/2j*bond.mod*bond_mod_D)*vecElm
                        ret_hamil[ind] += -1/4*bond_mod_J*vecElm
                        if(bond.mod == 1):
                            chachList[fripInd][1] += -2
                        else:
                            chachList[fripInd][1] += 2
                        chachList[fripInd][0] += 2
                        chachList[ind][0] += -1
                    elif(spin1 == 0 and spin2 == 1):
                        ret_hamil[fripInd] += (1/2*bond_mod_J+1/2j*bond.mod*bond_mod_D)*vecElm
                        ret_hamil[ind] += -1/4*bond_mod_J*vecElm
                        if(bond.mod == 1):
                            chachList[fripInd][1] += 2
                        else:
                            chachList[fripInd][1] += -2
                        chachList[fripInd][0] += 2
                        chachList[ind][0] += -1
                    else:
                        print("この状態はないよ")
            chach[ind] = chachList
    return ret_hamil

def main():
    global SITE_NUM,States,idToInd,chach
    #inputfileからの読み込み
    setter()
    #出力ファイルの準備
    logFile = open(f"./output/{inputFile}_J={bond_mod_J}_D={bond_mod_D}.txt","w")
    logFile.write("S_Z^2,E/j\n")
    #アップスピンとダウンスピンの数の差分DIFF_UP_DOWN_SPINに対して最低Eを取る(E(-1) = E(1))より、半分は割愛
    #DIFF_UP_DUWN_SPIN=SITE_NUMのときはすべてのスピンがそろっているので、基底になりえないから割愛
    programStart = datetime.datetime.now()
    for DIFF_UP_DOWN_SPIN in range(-SITE_NUM+2,1,2):
        startTime = datetime.datetime.now()
        #状態と逆変換を初期化
        States = []
        idToInd = {}
        chach = {}
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
        print("create states")
        hamilton_matrix = LinearOperator(
            shape=(ALL_STATE_NUM, ALL_STATE_NUM),
            matvec=Hamil,
            dtype=np.complex128
        )
        print("LinerOperator defaned")
        vals, vecs = eigsh(hamilton_matrix, k=1, which='SA') 
        print(f'MIN_ENEGY is {vals[0]}')
        print(datetime.datetime.now()-startTime)
        size_of_chach = sys.getsizeof(chach)
        print(f"chach配列のサイズ: {size_of_chach} バイト")
        logFile.write(f"{(DIFF_UP_DOWN_SPIN/2)**2},{vals[0]}\n")
    print(datetime.datetime.now()-programStart)
    logFile.close()

def statePrint():
    lis = []
    for s in States:
        binString = bin(s)[2:]
        lis.append("0"*(SITE_NUM - len(binString)) + binString)
    print(lis)
main()