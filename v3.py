import numpy as np
from functools import reduce
from itertools import combinations
import math
#クラス定義***************************
class Site():
    def __init__(self,Id,bond,direct):
        #Idはこのサイトの番号、bondはこれと結合しているサイトIdの配列(self.id<bond.id),directは方向(→を0とする度数法)
        if(len(bond) != len(direct)):
            raise ValueError("bond correspondence invaild")
        self.id = Id
        self.bond = bond
        self.direct = direct
    def __repr__(self):
        return f"State(id={self.id},bond = {self.bond})"
class State():
    def __init__(self,data):
        #dataは0,1配列でdata[i]はId=iのサイトの1,0を表す
        self.data = data
    def __eq__(self,other):
        return self.data == other.data
    def __repr__(self):
        return f"State(data={self.data})"
#定数***********************************
Sites = [
    Site(0,[1,15],[0,0]),
    Site(1,[2],[0]),
    Site(2,[3],[0]),
    Site(3,[4],[0]),
    Site(4,[5],[0]),
    Site(5,[6],[0]),
    Site(6,[7],[0]),
    Site(7,[8],[0]),
    Site(8,[9],[0]),
    Site(9,[10],[0]),
    Site(10,[11],[0]),
    Site(11,[12],[0]),
    Site(12,[13],[0]),
    Site(13,[14],[0]),
    Site(14,[15],[0]),
    Site(15,[],[])
]
SITE_NUM = len(Sites)
S_TOT = 0
ALL_STATE_NUM = math.comb(SITE_NUM,(S_TOT+SITE_NUM)//2)
bond_mod_J = 1
bond_mod_D = 1
#変数********************************
#Stateの配列
States = []
hamil = []
#関数定義******************************
def spinFrip(spin):
    if(spin == 1):
        return 0
    return 1
def Hamil(state):
    modified_state = [0]*ALL_STATE_NUM
    in_prod_ans = inProd(state)
    for [mod_state,mag] in in_prod_ans:
        modified_state[States.index(mod_state)] += mag*bond_mod_J
    out_prod_ans = outProd(state)
    for [mod_state,mag] in out_prod_ans:
        modified_state[States.index(mod_state)] += mag*bond_mod_D
    return modified_state

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
            copy_spins[site],copy_spins[bond_site] = copy_spins[bond_site],copy_spins[site]
            modified_state = State(copy_spins)
            ans.append([modified_state,1/4])
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
            ans.append([state,mag])
    return ans
#外積
def outProd(state):
    ans = []
    spins = state.data
    conv = {(1,1):0,(1,0):1j/4,(0,1):-1j/4,(0,0):0}
    for site in range(SITE_NUM):
        for bond_site in Sites[site].bond:
            copy_spins = spins.copy()
            copy_spins[site],copy_spins[bond_site] = copy_spins[bond_site],copy_spins[site]
            modified_state = State(copy_spins)
            ans.append([modified_state,conv[(spins[site],spins[bond_site])]])
    return ans

def main():
    print(f"this hamiltonian is {ALL_STATE_NUM} squred")
    for ups in combinations(range(SITE_NUM),(S_TOT+SITE_NUM)//2):
        spins = [0]*SITE_NUM
        for site in ups:
            spins[site] = 1
        States.append(State(spins))
    complete = 0
    disp = 1
    for state in States:
        mod_list = Hamil(state)
        hamil.append(mod_list)
        complete += 1
        if(complete/ALL_STATE_NUM >= disp/100):
            print(f"hamiltonian {disp}% complate")
            disp += 1
        #print(state,mod_list)
    #print(hamil)
    print("hamiltonial was created")
    [koyuuti,vec] = np.linalg.eig(hamil)
    koyuuti = sorted(koyuuti,key = lambda i:i.real)
    for i in range(len(koyuuti)):
        if(koyuuti[i].imag > 1):
            print("imag part exest")
    dis = 0
    print(len(koyuuti)," exest")
    while len(koyuuti)>dis:
        print(dis,koyuuti[dis].real)
        dis += 10
    print(len(koyuuti)-1,koyuuti[-1].real)
main()