import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import eigs
from scipy.linalg import eig
from functools import reduce
from itertools import combinations
import math
import datetime

#クラス定義***************************
class Site():
    def __init__(self,Id,bond,mod):
        #Idはこのサイトの番号、bondはこれと結合しているサイトIdの配列(self.id<bond.id),modは定数(+上向き)
        if(len(bond) != len(mod)):
            raise ValueError("bond correspondence invaild")
        self.id = Id
        self.bond = bond
        self.mod = mod
    def __repr__(self):
        return f"State(id={self.id},bond = {self.bond})"
class State():
    def __init__(self,data):
        #dataは0,1配列でdata[i]はId=iのサイトの1,0を表す
        self.data = data
        Id = 0
        for spin in self.data:
            Id = Id*2+spin
        self.id = Id
    def __eq__(self,other):
        return self.data == other.data
    def __repr__(self):
        return f"State(data={self.data})"
#定数***********************************
bond_mod_J = 1
bond_mod_D = 0.1
accuracy = 1e-2 #精度
Sites = []
SITE_NUM = 0
inputFile = "N=21"
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
    stateFile = open("./input/" + inputFile+".txt","r")
    num = int(stateFile.readline())
    for line in range(num):
        print(line)
        data = stateFile.readline().replace("\n", "")
        if(data != ""):
            bond = list(map(int,data.split(",")))
        else:
            bond = []
        data = stateFile.readline().replace("\n", "")
        if(data != ""):
            mod = list(map(int,data.split(",")))
        else:
            mod = []
        Sites.append(Site(line,bond,list(map(lambda moddata: moddata*bond_mod_D,mod))))
    stateFile.close()
#ハミルトニアン操作
def Hamil(state):
    in_prod_ans = inProd(state)
    out_prod_ans = outProd(state)
    return in_prod_ans + out_prod_ans
#spin=>Id
def spinToId(spins):
    Id = 0
    for spin in spins:
        Id = Id*2+spin
    return Id
#内積
def inProd(state):
    #stateはState型
    #この配列には[State(),係数]を入れる
    ans = []
    x_ans = xyprod(state)
    ans += x_ans
    z_ans = zprod(state)
    ans += z_ans
    return ans
#x,yの内積
def xyprod(state):
    ans = []
    spins = state.data
    for site in range(SITE_NUM):
        for bond_site in Sites[site].bond:
            copy_spins = spins.copy()
            if(copy_spins[site] != copy_spins[bond_site]):
                copy_spins[site],copy_spins[bond_site] = copy_spins[bond_site],copy_spins[site]
                modified_state = States[idToInd[spinToId(copy_spins)]]
                ans.append([modified_state,1/2*bond_mod_J])
    return ans
#zの内積
def zprod(state):
    ans = []
    spins = state.data
    for site in range(SITE_NUM):
        for bond_site in Sites[site].bond:
            mag = 1
            if(spins[site] == 1):
                mag *= 1/2
            else:
                mag *= -1/2
            if(spins[bond_site] == 1):
                mag *= 1/2
            else:
                mag *= -1/2
            ans.append([state,mag*bond_mod_J])
    return ans
#外積
def outProd(state):
    ans = []
    spins = state.data
    conv = {(1,0):1j/2,(0,1):-1j/2}
    for site in range(SITE_NUM):
        for bond_ind in range(len(Sites[site].bond)):
            copy_spins = spins.copy()
            if(copy_spins[site] != copy_spins[Sites[site].bond[bond_ind]]):
                copy_spins[site],copy_spins[Sites[site].bond[bond_ind]] = copy_spins[Sites[site].bond[bond_ind]],copy_spins[site]
                modified_state = States[idToInd[spinToId(copy_spins)]]
                ans.append([modified_state,conv[(spins[site],spins[Sites[site].bond[bond_ind]])]*Sites[site].mod[bond_ind]])
    return ans

def main():
    startTime = datetime.datetime.now()
    global SITE_NUM
    #inputfileからの読み込み
    setter()
    SITE_NUM = len(Sites)
    #出力ファイルの準備
    logFile = open(f"./output/N={SITE_NUM}_J={bond_mod_J}_D={bond_mod_D}.txt","w")
    logFile.write("S_Z^2,E/j\n")
    min_enegy = float("inf")
    #全S_TOTに対して最低Eを取る(E(-1) = E(1))より、半分は割愛
    for S_TOT in range(SITE_NUM,-1,-2):
        #状態と逆変換を初期化
        States.clear()
        idToInd.clear()
        print(f'S_TOT is {S_TOT}')
        ALL_STATE_NUM = math.comb(SITE_NUM,(S_TOT+SITE_NUM)//2)
        print(f"this hamiltonian is {ALL_STATE_NUM} squred")
        StateInd = 0
        #状態と逆変換を列挙
        for ups in combinations(range(SITE_NUM),(S_TOT+SITE_NUM)//2):
            spins = [0]*SITE_NUM
            for site in ups:
                spins[site] = 1
            States.append(State(spins))
            idToInd[State(spins).id] = StateInd
            StateInd += 1
        #どこまで終わったか出力するため
        complete = 0
        disp = 10
        dataDict = {}
        #ハミルトニアン形成{(元のId,変換後のID) = 係数}
        for state in States:
            col_data = Hamil(state)
            rowId = idToInd[state.id]
            for _col,_data in col_data:
                colId = idToInd[_col.id]
                if((rowId,colId) in dataDict):
                    dataDict[(rowId,colId)] += _data
                else:
                    dataDict[(rowId,colId)] = _data
            complete += 1
            if(complete/ALL_STATE_NUM >= disp/100):
                print(f"hamiltonian {math.floor(complete*10//ALL_STATE_NUM)*10}% complate")
                disp = math.floor(complete*10//ALL_STATE_NUM)*10 + 10
            #print(state,mod_list)
        print("convart was finished")
        #疎行列用にコンバート
        max_nnz = len(dataDict)
        rows = np.empty(max_nnz, dtype=int)
        cols = np.empty(max_nnz, dtype=int)
        datas = np.empty(max_nnz, dtype=complex)
        ptr = 0
        for key,data in dataDict.items():
            rows[ptr] = key[0]
            cols[ptr] = key[1]
            datas[ptr] = data
            ptr += 1
        print("create list data")
        #結局ハミルトニアンが小さいなら蜜行列を使った方が早いので実装
        if(ALL_STATE_NUM <= 100):
            hamil = coo_matrix((datas, (rows, cols)), shape=(ALL_STATE_NUM, ALL_STATE_NUM)).tocsr()
            print("hamiltonial was created")
            eigvals, eigvecs = eig(hamil.toarray())
            mineigval = float("inf")
            for enegy in eigvals:
                mineigval = min(mineigval,enegy.real)
            min_enegy = min(min_enegy,mineigval) 
            print(S_TOT,"環境下での最小実部の固有値:", math.floor(mineigval/accuracy*10)*accuracy/10)
        #疎行列の時はこっち
        else:
            hamil = coo_matrix((datas, (rows, cols)), shape=(ALL_STATE_NUM, ALL_STATE_NUM)).tocsr()
            print("hamiltonial was created")
            eigvals, eigvecs = eigs(hamil, k=1, which='SR',tol=accuracy)
            mineigval = float("inf")
            for enegy in eigvals:
                mineigval = min(mineigval,enegy.real)
            min_enegy = min(min_enegy,mineigval) 
            print(S_TOT,"環境下での最小実部の固有値:", math.floor(mineigval/accuracy*10)*accuracy/10)
        logFile.write(f"{(S_TOT/2)**2},{eigvals[0].real}\n")
    print("全ての状態での最小エネルギー:",min_enegy)
    logFile.close()
    print(datetime.datetime.now()-startTime)

main()