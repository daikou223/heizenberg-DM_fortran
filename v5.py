import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import eigs
from scipy.linalg import eig
from functools import reduce
from itertools import combinations
import math
import datetime
from bekizyou import koyuuti

#クラス定義***************************
class Site():
    def __init__(self,Id,bond,mod):
        #Idはこのサイトの番号、bondはこれと結合しているサイトIdの配列(self.id<bond.id),modは定数(+上向き)
        self.id = Id
        self.bond = bond
        self.mod = mod
    def __repr__(self):
        return f"Site(id={self.id},bond = {self.bond},mod = {self.mod})"
#定数***********************************
bond_mod_J = 1
bond_mod_D = 0.1
accuracy = 1e-2 #精度
Sites = []
SITE_NUM = 0
inputFile = "N=16_test"
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
    return num
def get_bit(num,i):
    return num >> (SITE_NUM-1-i) & 1
#ハミルトニアン操作
def Hamil(state):
    in_prod_ans = inProd(state)
    out_prod_ans = outProd(state)
    return in_prod_ans + out_prod_ans
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
    spins = state
    for site in range(SITE_NUM):
        for bond_site in Sites[site].bond:
            copy_spins = spins
            if((get_bit(copy_spins,site) & 1) != (get_bit(copy_spins,bond_site) & 1)):
                copy_spins ^= 1<<(SITE_NUM-1-site) | 1<<(SITE_NUM-1-bond_site)
                ans.append([copy_spins,1/2*bond_mod_J])
    return ans
#zの内積
def zprod(state):
    ans = []
    spins = state
    for site in range(SITE_NUM):
        for bond_site in Sites[site].bond:
            mag = 1
            if(get_bit(spins,site) == 1):
                mag *= 1/2
            else:
                mag *= -1/2
            if(get_bit(spins,bond_site) == 1):
                mag *= 1/2
            else:
                mag *= -1/2
            ans.append([spins,mag*bond_mod_J])
    return ans
#外積
def outProd(state):
    ans = []
    spins = state
    conv = {(1,0):-1j/2,(0,1):+1j/2}
    for site in range(SITE_NUM):
        for bond_ind,bond in enumerate(Sites[site].bond):
            copy_spins = spins
            if(get_bit(copy_spins,site) != get_bit(copy_spins,bond)):
                copy_spins ^= 1 << (SITE_NUM-1-site) | 1 << (SITE_NUM-1-bond)
                ans.append([copy_spins,conv[(get_bit(spins,site),get_bit(spins,bond))]*Sites[site].mod[bond_ind]])
    return ans

def main():
    startTime = datetime.datetime.now()
    global SITE_NUM,States,idToInd
    #inputfileからの読み込み
    num = setter()
    SITE_NUM = num
    #出力ファイルの準備
    logFile = open(f"./output/{inputFile}_J={bond_mod_J}_D={bond_mod_D}.txt","w")
    logFile.write("S_Z^2,E/j\n")
    min_enegy = float("inf")
    #全S_TOTに対して最低Eを取る(E(-1) = E(1))より、半分は割愛
    for S_TOT in range(-SITE_NUM,1,2):
        #状態と逆変換を初期化
        States = []
        idToInd = {}
        print(f'S_TOT is {S_TOT}')
        ALL_STATE_NUM = math.comb(SITE_NUM,(S_TOT+SITE_NUM)//2)
        print(f"this hamiltonian is {ALL_STATE_NUM} squred")
        StateInd = 0
        #状態と逆変換を列挙
        for ups in combinations(range(SITE_NUM),(S_TOT+SITE_NUM)//2):
            spins = 0
            for i in ups:
                spins |= (1<<(SITE_NUM-1-i))
            States.append(spins)
            idToInd[spins] = StateInd
            StateInd += 1
        #どこまで終わったか出力するため
        complete = 0
        disp = 10
        dataDict = []
        dataNum = []
        #ハミルトニアン形成{(元のId,変換後のID) = 係数}
        for state in States:
            col_data = Hamil(state)
            rowId = idToInd[state]
            dataDict.append({})
            dataNum.append(0)
            for _col,_data in col_data:
                colId = idToInd[_col]
                if(colId in dataDict[rowId]):
                    dataDict[rowId][colId] += _data
                else:
                    dataDict[rowId][colId] = _data
                    dataNum[-1] += 1
            #進捗表示
            complete += 1
            if(complete/ALL_STATE_NUM >= disp/100):
                print(f"hamiltonian {math.floor(complete*10//ALL_STATE_NUM)*10}% complate")
                disp = math.floor(complete*10//ALL_STATE_NUM)*10 + 10
        #疎行列用にコンバート
        _rows = np.arange(ALL_STATE_NUM)
        rows = np.repeat(_rows, dataNum)
        total_length = sum(dataNum)
        cols = np.empty(total_length, dtype=int)
        datas = np.empty(total_length, dtype=complex)
        ptr = 0
        rep = 0
        for wrapper in dataDict:
            # 一度にまとめて取得
            col_vals, data_vals = zip(*wrapper.items())
            cols[ptr:ptr+dataNum[rep]] = col_vals
            datas[ptr:ptr+dataNum[rep]] = data_vals
            ptr += dataNum[rep]
            rep += 1
        #結局ハミルトニアンが小さいなら蜜行列を使った方が早いので実装
        #if(ALL_STATE_NUM <= 100):
        hamil = coo_matrix((datas, (rows, cols)), shape=(ALL_STATE_NUM, ALL_STATE_NUM)).tocsr()
        print("hamiltonial was created")
        eigval = koyuuti(hamil,ALL_STATE_NUM)
        # eigvals, eigvecs = eig(hamil.toarray())
        # mineigval = float("inf")
        # for enegy in eigvals:
        #     mineigval = min(mineigval,enegy.real)
        # min_enegy = min(min_enegy,mineigval) 
        print(S_TOT,"環境下での最小実部の固有値:", math.floor(eigval/accuracy*10)*accuracy/10)
        #疎行列の時はこっち
        # else:
        #     hamil = coo_matrix((datas, (rows, cols)), shape=(ALL_STATE_NUM, ALL_STATE_NUM)).tocsr()
        #     print("hamiltonial was created")
        #     eigvals, eigvecs = eigs(hamil, k=1, which='SR',tol=accuracy)
        min_enegy = min(min_enegy,eigval.real) 
        #     print(S_TOT,"環境下での最小実部の固有値:", math.floor(eigvals[0].real/accuracy*10)*accuracy/10)
        logFile.write(f"{(S_TOT/2)**2},{eigval.real}\n")
    print("全ての状態での最小エネルギー:",min_enegy)
    logFile.close()
    print(datetime.datetime.now()-startTime)

main()