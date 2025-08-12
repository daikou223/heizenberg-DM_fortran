import numpy as np
from itertools import combinations
import math
import datetime
from scipy.sparse.linalg import lobpcg,LinearOperator, eigsh, gmres
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
inputFile = f"N={SITE_NUM}"
#変数********************************
#Stateの配列
States = []
idToInd = {}
hamil = []
rows = []
cols = []
datas = []
UP_NUM = 0
beforeData = np.array([])
beforeEnegy = float("inf")
enegyData = [float("inf"),float("inf"),float("inf")]
predictEnegy = float("inf")
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
    global SITE_NUM,UP_NUM,beforeData
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

def realyOper(vec):
    return inplod(vec,Hamil(vec))

def nableRealy(vec,r):
    return (Hamil(vec)-min(r,beforeEnegy,predictEnegy)*vec)

def inplod(vec1,vec2):
    inplod_ = np.vdot(vec1, vec2)
    return inplod_

def main():
    global SITE_NUM,States,idToInd,UP_NUM, beforeData, beforeEnegy, enegyData, predictEnegy
    #inputfileからの読み込み
    setter()
    #出力ファイルの準備
    logFile = open(f"./output/{os.path.basename(__file__)}_{inputFile}_J={bond_mod_J}_D={bond_mod_D}.txt","w")
    logFile.write("S_Z^2,E/j\n")
    programStart = datetime.datetime.now()
    #アップスピンとダウンスピンの数の差分DIFF_UP_DOWN_SPINに対して最低Eを取る(E(-1) = E(1))より、半分は割愛
    #DIFF_UP_DUWN_SPIN=SITE_NUMのときはすべてのスピンがそろっているので、基底になりえないから割愛
    for UP_NUM in range(1,SITE_NUM,1):
        DOWN_NUM = SITE_NUM-UP_NUM
        #逆のときは対称性より計算しなくてよい
        if(UP_NUM<=DOWN_NUM):
            enegyData = [float("inf"),float("inf"),float("inf")]
            startTime = datetime.datetime.now()
            #状態と逆変換を初期化
            beforeIdToIndex = idToInd
            States = []
            idToInd = {}
            print(f'UP_NUM is {UP_NUM}')
            ALL_STATE_NUM = math.comb(SITE_NUM,UP_NUM)
            print(f"this hamiltonian is {ALL_STATE_NUM} squred")
            StateInd = 0
            X = np.random.normal(scale = 1e-3,size = ALL_STATE_NUM)
            #状態と逆変換を列挙
            for ups in combinations(range(SITE_NUM),UP_NUM):
                spins = 0
                alikeSpins = 0
                for i in ups:
                    spins |= (1<<i)
                for i in range(len(ups)-1):
                    alikeSpins |= 1 << ups[i]
                States.append(spins)
                idToInd[spins] = StateInd
                if(alikeSpins in beforeIdToIndex):
                    X[StateInd] = beforeData[beforeIdToIndex[alikeSpins]]
                StateInd += 1
            X = X/inplod(X,X)**0.5
            Attenuation = 1
            beforeR = float("inf")
            r = 0
            loopcount = 0
            dispind = 0
            while 1:
                beforeR = r
                r = realyOper(X)
                hamilX = Hamil(X)
                nablaR = nableRealy(X,r)
                hamilNablaR = Hamil(nablaR)
                #二部探索を行う
                maxAttenuation = Attenuation*2
                minAttenuation = 0
                nablaRHamilX = inplod(nablaR,hamilX)
                xHamilNablaR = inplod(X,hamilNablaR)
                nablaRHamilNablaR = inplod(nablaR,hamilNablaR)
                nablaRX = inplod(nablaR,X)
                nablaRNablaR = inplod(nablaR,nablaR)
                while maxAttenuation-minAttenuation > 10**-3:
                    Attenuation = (maxAttenuation+minAttenuation)/2
                    DAttenuation = (maxAttenuation+minAttenuation)/2+10**-8
                    newR = (r + Attenuation*(-nablaRHamilX - xHamilNablaR) + Attenuation**2*nablaRHamilNablaR)/(1-2*Attenuation*nablaRX+Attenuation**2*nablaRNablaR)
                    DnewR = (r + DAttenuation*(-nablaRHamilX - xHamilNablaR) + DAttenuation**2*nablaRHamilNablaR)/(1-2*DAttenuation*nablaRX+DAttenuation**2*nablaRNablaR)
                    gradient = (DnewR - newR) / (DAttenuation - Attenuation)
                    if(gradient.real > 0):
                        maxAttenuation = Attenuation
                    else:
                        minAttenuation = Attenuation
                X = (X-Attenuation*nablaR)/(1-2*Attenuation*nablaRX+Attenuation**2*nablaRNablaR)**0.5
                loopcount += 1
                if(loopcount == 2**dispind):
                    print(Attenuation)
                    print((beforeR-r).real,r.real,inplod(nablaR,nablaR).real,dispind)
                    enegyData[2] = enegyData[1]
                    enegyData[1] = enegyData[0]
                    enegyData[0] = r.real
                    if(enegyData[2] != float("inf")):
                        lamN = enegyData[2]
                        lam2N = enegyData[1]
                        lam4N = enegyData[0]
                        k = (lamN-lam2N)/(lam2N-lam4N)
                        beta = (-1+(1+4/k)**0.5)/2
                        predictEnegy = lamN-((lamN-lam2N)**2*beta)/(2*beta*(lamN-lam2N)-(lam2N-lam4N))
                        if(abs(predictEnegy-r.real)<10**-2):
                            break
                    dispind += 1
                if(inplod(nablaR,nablaR).real < 0.001):
                    break
            print(f"MINENEGY = {r.real}")
            beforeData = X
            beforeEnegy = r.real
            logFile.write(f"{(UP_NUM*1/2-DOWN_NUM*1/2)**2},{r.real},{str(datetime.datetime.now()-startTime)}\n")
    logFile.write(str(datetime.datetime.now()-programStart))
    logFile.close()

main()