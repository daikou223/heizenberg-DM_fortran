def getIncludeSite(Apos,Bpos,target_site = float("inf")):
    sitePos = set()
    minuteApos = [Apos[0]*2,Apos[1]*2]
    minuteBpos = [Bpos[0]*2,Bpos[1]*2]
    #平行の場合は即座に返却
    if(minuteApos[1]*minuteBpos[0]==minuteApos[0]*minuteBpos[1]):
        return 0
    startLine = min(0,minuteApos[0],minuteBpos[0],minuteApos[0]+minuteBpos[0])
    endLine = max(0,minuteApos[0],minuteBpos[0],minuteApos[0]+minuteBpos[0])
    startRow = min(0,minuteApos[1],minuteBpos[1],minuteApos[1]+minuteBpos[1])
    endRow = max(0,minuteApos[1],minuteBpos[1],minuteApos[1]+minuteBpos[1])
    #これらは必ず偶数
    count = 0
    for aline in range(startLine,endLine+1):
        for arow in range(startRow,endRow):
            isIncludeFrag = isInclude(aline,arow,minuteApos,minuteBpos)
            if(isIncludeFrag == 1):
                count += 1
                sitePos.add((aline,arow))
                if(count > target_site):
                    return target_site+1,set()
    return count,sitePos

def isInclude(aline,arow,minuteApos,minuteBpos):
    #0:そんなサイトがない 1:内部に存在 2:
    if((aline%2 == 0 and arow%2 == 0) or (arow%4 == 1 and aline%2 == 1)):
        l = (minuteApos[1]*aline-minuteApos[0]*arow)
        if(l < 0):
            signl = -1
        else:
            signl = 1
        k = (minuteBpos[1]*aline-minuteBpos[0]*arow)
        if(k < 0):
            signk = -1
        else:
            signk = 1
        if((l == 0 and 0 <= k*signk <= (minuteApos[0]*minuteBpos[1]-minuteApos[1]*minuteBpos[0])*signk) or (k == 0 and 0 <= l*signl <= (minuteApos[1]*minuteBpos[0]-minuteApos[0]*minuteBpos[1])*signl) or 0 <= l*signl < (minuteApos[1]*minuteBpos[0]-minuteApos[0]*minuteBpos[1])*signl and 0 <= k*signk < (minuteApos[0]*minuteBpos[1]-minuteApos[1]*minuteBpos[0])*signk):
            return 1
        elif(l*signl >= (minuteApos[1]*minuteBpos[0]-minuteApos[0]*minuteBpos[1])*signl):
            if(k/(minuteApos[0]*minuteBpos[1]-minuteApos[1]*minuteBpos[0]) < 0):
                return 2
            elif(k*signk < (minuteApos[0]*minuteBpos[1]-minuteApos[1]*minuteBpos[0])*signk):
                return 4
            else:
                return 3
        elif(l/(minuteApos[1]*minuteBpos[0]-minuteApos[0]*minuteBpos[1]) <= 0):
            if(k/(minuteApos[0]*minuteBpos[1]-minuteApos[1]*minuteBpos[0]) < 0):
                return 8
            elif(k*signk < (minuteApos[0]*minuteBpos[1]-minuteApos[1]*minuteBpos[0])*signk):
                return 7
            else:
                return 6
        else:
            if(k/(minuteApos[0]*minuteBpos[1]-minuteApos[1]*minuteBpos[0]) < 0):
                return 9
            elif(k*signk < (minuteApos[0]*minuteBpos[1]-minuteApos[1]*minuteBpos[0])*signk):
                return 5
    return 0

def main():
    target_site = 21
    upperrange = 5 
    #広い範囲で大変すぎたのでレンジをいい感じに設定
    for aline in range(-upperrange,upperrange+1):
        print(aline)
        for bline in range(-upperrange,upperrange+1):
            for arow in range(-upperrange,upperrange+1):
                for brow in range(arow,upperrange+1):
                    #きれいにトーラスになる条件
                    if(aline%2 == arow%2 and bline%2 == brow%2):
                        includeSite,sitePos = getIncludeSite([aline,arow],[bline,brow],target_site)
                        if(includeSite == target_site):
                            posToId = {}
                            index = 0
                            inputFile = open(f'./../input/${(aline,arow)},${(bline,brow)}_N=${target_site}.txt',"w")
                            for sitePo in sitePos:
                                posToId[sitePo] = index
                                index += 1
                            for sitePo in sitePos:
                                if(aline % 2 == 0):
                                    newLine = sitePo[0]-2
                                    newRow = sitePo[1]
                                    BSite = posToId[sitePo]
                                    if((newLine,newRow) in posToId):
                                        ASite = posToId[(newLine,newRow)]
                                        inputFile.write(f'f{BSite} {ASite} 1')
                                        continue
                            #上の処理よくわかんなくなってきたからあとよろしく~   
                            print((aline,arow),(bline,brow))

main()