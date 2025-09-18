from fractions import Fraction
import os

def getIncludeSite(Apos,Bpos,target_site = float("inf")):
    sitePos = set()
    #平行の場合は即座に返却
    if(Apos[1]*Bpos[0]==Apos[0]*Bpos[1]):
        return 0,set()
    startLine = min(0,Apos[0],Bpos[0],Apos[0]+Bpos[0])
    endLine = max(0,Apos[0],Bpos[0],Apos[0]+Bpos[0])
    startRow = min(0,Apos[1],Bpos[1],Apos[1]+Bpos[1])
    endRow = max(0,Apos[1],Bpos[1],Apos[1]+Bpos[1])
    #これらは必ず偶数
    count = 0
    for aline in range(startLine,endLine+1):
        for arow in range(startRow,endRow+1):
            isSiteFlag = False
            if(aline %2 == 0 and arow % 2 == 0):
                isSiteFlag = True
            if(aline %4 == 1 and arow % 4 == 1):
                isSiteFlag = True
            if(aline %4 == 3 and arow % 4 == 3):
                isSiteFlag = True
            if(isSiteFlag):
                k,l = decomposeVerter(aline,arow,Apos,Bpos)
                if((0 <= l < 1 and 0 <= k < 1)):
                    count += 1
                    sitePos.add((aline,arow))
                    if(count > target_site):
                        return target_site+1,set()
    return count,sitePos

def decomposeVerter(aline,arow,Apos,Bpos):
    if(Apos[1]*Bpos[0]-Apos[0]*Bpos[1] != 0 and Apos[0]*Bpos[1]-Apos[1]*Bpos[0]):
        l = Fraction((Apos[1]*aline-Apos[0]*arow),(Apos[1]*Bpos[0]-Apos[0]*Bpos[1]))
        k = Fraction((Bpos[1]*aline-Bpos[0]*arow),(Apos[0]*Bpos[1]-Apos[1]*Bpos[0]))
        return k,l
    return float("inf"),float("inf")

class space():
    def __init__(self,aline,arow,bline,brow):
        self.aline = aline
        self.arow = arow
        self.bline = bline
        self.brow = brow
        #ベクトルの長さの和が小さい方を優先して先に計算
        point = aline**2+arow**2+bline**2+brow**2
        if(aline < 0):
            point *= 1.1
        if(arow < 0):
            point *= 1.1
        if(bline < 0):
            point *= 1.1
        if(brow < 0):
            point *= 1.1
        self.point = point

def main():
    target_site = 6
    upperrange = 6
    passedList = []
    os.makedirs(f'./input/N{target_site}', exist_ok=True)
    #広い範囲で大変すぎたのでレンジをいい感じに設定
    for aline in range(-2*upperrange,2*upperrange+1,2):
        for bline in range(-2*upperrange,2*upperrange+1,2):
            for arow in range(-2*upperrange,2*upperrange+1,2):
                for brow in range(arow,2*upperrange+1,2):
                    #きれいにトーラスになる条件
                    if(aline%4 == arow % 4 and bline%4 == brow%4):
                        includeSite,sitePos = getIncludeSite([aline,arow],[bline,brow],target_site)
                        if(includeSite == target_site):
                            posToId = {}
                            #各座標をナンバリング
                            index = 0
                            inputFile = open((f'./input/N{target_site}/{(aline,arow)}_{(bline,brow)}_N{target_site}.txt').replace(" ",""),"w")
                            for sitePo in sitePos:
                                posToId[sitePo] = index
                                index += 1
                            #ボンド結合がある場合
                            for sitePo in sitePos:
                                ASiteId = posToId[sitePo]
                                if(sitePo[0] %4 == 0 and sitePo[1] % 4 == 0):
                                    newLine = sitePo[0]
                                    newRow = sitePo[1]-2
                                    k,l = decomposeVerter(newLine,newRow,[aline,arow],[bline,brow])
                                    l = l % 1
                                    k = k % 1
                                    newLine = k*aline + l*bline
                                    newRow = k*arow + l*brow
                                    BSiteId = posToId[(newLine,newRow)]
                                    inputFile.write(f'{ASiteId} {BSiteId} -1\n')

                                    newLine = sitePo[0]+1
                                    newRow = sitePo[1]+1
                                    k,l = decomposeVerter(newLine,newRow,[aline,arow],[bline,brow])
                                    l = l % 1
                                    k = k % 1
                                    newLine = k*aline + l*bline
                                    newRow = k*arow + l*brow
                                    BSiteId = posToId[(newLine,newRow)]
                                    inputFile.write(f'{ASiteId} {BSiteId} 1\n')
                                if(sitePo[0] %4 == 0 and sitePo[1] % 4 == 2):
                                    newLine = sitePo[0]
                                    newRow = sitePo[1]-2
                                    k,l = decomposeVerter(newLine,newRow,[aline,arow],[bline,brow])
                                    l = l % 1
                                    k = k % 1
                                    newLine = k*aline + l*bline
                                    newRow = k*arow + l*brow
                                    BSiteId = posToId[(newLine,newRow)]
                                    inputFile.write(f'{ASiteId} {BSiteId} 1\n')

                                    newLine = sitePo[0]-1
                                    newRow = sitePo[1]+1
                                    k,l = decomposeVerter(newLine,newRow,[aline,arow],[bline,brow])
                                    l = l % 1
                                    k = k % 1
                                    newLine = k*aline + l*bline
                                    newRow = k*arow + l*brow
                                    BSiteId = posToId[(newLine,newRow)]
                                    inputFile.write(f'{ASiteId} {BSiteId} -1\n')
                                if(sitePo[0] %4 == 2 and sitePo[1] % 4 == 0):
                                    newLine = sitePo[0]
                                    newRow = sitePo[1]-2
                                    k,l = decomposeVerter(newLine,newRow,[aline,arow],[bline,brow])
                                    l = l % 1
                                    k = k % 1
                                    newLine = k*aline + l*bline
                                    newRow = k*arow + l*brow
                                    BSiteId = posToId[(newLine,newRow)]
                                    inputFile.write(f'{ASiteId} {BSiteId} 1\n')

                                    newLine = sitePo[0]-1
                                    newRow = sitePo[1]+1
                                    k,l = decomposeVerter(newLine,newRow,[aline,arow],[bline,brow])
                                    l = l % 1
                                    k = k % 1
                                    newLine = k*aline + l*bline
                                    newRow = k*arow + l*brow
                                    BSiteId = posToId[(newLine,newRow)]
                                    inputFile.write(f'{ASiteId} {BSiteId} -1\n')
                                if(sitePo[0] %4 == 2 and sitePo[1] % 4 == 2):
                                    newLine = sitePo[0]
                                    newRow = sitePo[1]-2
                                    k,l = decomposeVerter(newLine,newRow,[aline,arow],[bline,brow])
                                    l = l % 1
                                    k = k % 1
                                    newLine = k*aline + l*bline
                                    newRow = k*arow + l*brow
                                    BSiteId = posToId[(newLine,newRow)]
                                    inputFile.write(f'{ASiteId} {BSiteId} -1\n')

                                    newLine = sitePo[0]+1
                                    newRow = sitePo[1]+1
                                    k,l = decomposeVerter(newLine,newRow,[aline,arow],[bline,brow])
                                    l = l % 1
                                    k = k % 1
                                    newLine = k*aline + l*bline
                                    newRow = k*arow + l*brow
                                    BSiteId = posToId[(newLine,newRow)]
                                    inputFile.write(f'{ASiteId} {BSiteId} 1\n')
                                if(sitePo[0] %2 == 1 and sitePo[1] % 2 == 1):
                                    newLine = sitePo[0]+1
                                    newRow = sitePo[1]+1
                                    k,l = decomposeVerter(newLine,newRow,[aline,arow],[bline,brow])
                                    l = l % 1
                                    k = k % 1
                                    newLine = k*aline + l*bline
                                    newRow = k*arow + l*brow
                                    BSiteId = posToId[(newLine,newRow)]
                                    inputFile.write(f'{ASiteId} {BSiteId} -1\n')

                                    newLine = sitePo[0]-1
                                    newRow = sitePo[1]+1
                                    k,l = decomposeVerter(newLine,newRow,[aline,arow],[bline,brow])
                                    l = l % 1
                                    k = k % 1
                                    newLine = k*aline + l*bline
                                    newRow = k*arow + l*brow
                                    BSiteId = posToId[(newLine,newRow)]
                                    inputFile.write(f'{ASiteId} {BSiteId} 1\n')
                            if(aline == 4 and arow == 0 and bline == 0 and brow == 4):
                                print(posToId)
                            inputFile.close()
                            passedList.append(space(aline,arow,bline,brow))
    if(len(passedList) > 0):
        ListFile = open(f'./input/N{target_site}/List.txt',"w")
        passedList.sort(key=lambda x:x.point)
        for i in passedList:
            ListFile.write(f'{i.aline},{i.arow},{i.bline},{i.brow}\n')
        ListFile.close()
    else:
        print("存在しません")
        os.rmdir(f'./input/N{target_site}')
main()
