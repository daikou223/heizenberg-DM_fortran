program baseEnegy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    !可変部
    real,parameter::J =1.,D = 0.1,accuracy = 1e-5
    integer,parameter::N = 18,maxSearch = 1
    integer,parameter::bond_Num = 2*n
    complex,parameter::reaser = 30
    real,parameter::debugNum = 0.01

    !非可変部
    character(len = 100)::inputPath,outputPath,chashPath,statePath,vectorPath
    character(len = 100)::tmpChar,tmpRead1,tmpRead2,tmpRead,ReadChash
    character(len = 10)::vectors(maxSearch,2,2)
    integer,allocatable::states(:)
    integer::Bonds(bond_Num,3)
    integer::ALL_STATE_NUM,valInd
    real::enegys(N),val,vals(N)
    integer:: v1,t1,t2,t3,t4,flag,restStatesNum
    integer::l1,l2,UP_NUM,DOWN_NUM
    integer, allocatable :: combList(:),restStates(:)
    complex, allocatable::newState(:),initState(:)
    integer::testnum
    complex,allocatable::ret_hamil(:),initData(:),testVec(:),testInitVec(:)
    integer :: start_count, end_count, rate, count_max, times(N),idToIndChash(0:2**N-1)
    logical::fileExists
    !ファイル書き込み用
    character(len=50) :: strN, strJ, strD, strU

    write(tmpChar,'(I2)') N
    inputPath = "./input/N" // trim(adjustl(tmpchar)) // "/List.txt"
    open(unit=10, file=inputPath, status="old")
        do l1 = 1,maxSearch
            read(10,'(A)') tmpRead1
            read(tmpRead1,*) t1,t2,t3,t4
            write(vectors(l1,1,1),'(I5)') t1
            write(vectors(l1,1,2),'(I5)') t2
            write(vectors(l1,2,1),'(I5)') t3
            write(vectors(l1,2,2),'(I5)') t4
        end do
    close(10)
    do vecInd = 1,maxSearch
        call system_clock(start_count, rate, count_max)
        vectorPath = "(" // trim(adjustl(vectors(vecInd,1,1))) // &
                     "," // trim(adjustl(vectors(vecInd,1,2))) // &
                     ")_(" // trim(adjustl(vectors(vecInd,2,1))) // &
                     "," // trim(adjustl(vectors(vecInd,2,2))) // &
                     ")" // "_N" // trim(adjustl(tmpchar)) // ".txt"
        statePath = "./input/N" // trim(adjustl(tmpchar)) // &
                     "/" // vectorPath 
        open(unit=10, file=statePath, status="old")
            do l2 = 1,bond_Num
                read(10,'(A)') tmpRead1
                read(tmpRead1,*) t1,t2,t3
                Bonds(l2,1) = t1
                Bonds(l2,2) = t2
                Bonds(l2,3) = t3
            end do
        close(10)
        do UP_NUM = 1,N
            print *,"UP_NUM = ",UP_NUM
            DOWN_NUM = N-UP_NUM
            !逆のときは対称性より計算しない
            if(UP_NUM <= DOWN_NUM) then
                ALL_STATE_NUM = nCr(N,UP_NUM)
                print *,"HamilSize ",ALL_STATE_NUM
                states = combinations(N,UP_NUM)
                restStates = combinations(N-2,UP_NUM-1)
                restStatesNum = nCr(N-2,UP_NUM-1)
                do L2 = 1,ALL_STATE_NUM
                    idToIndChash(states(L2)) = idToInd(states(L2),N,UP_NUM)
                end do
                print *,"create Chash"
                initState = random_vec_gene(ALL_STATE_NUM)
                print *,"start raily method"
                val = raily(initState,ALL_STATE_NUM,accuracy,Bonds,J,D,UP_NUM,restStates,restStatesNum,chashPath,idToIndChash)
                print *,val
                vals(UP_NUM) = val
                valInd = UP_NUM
                call system_clock(end_count)
                times(UP_NUM) = real(end_count - start_count) / real(rate)
                beforeCount = end_count
            end if
        end do

        call system_clock(end_count)
        print *, "Execution time (wall): ", real(end_count - start_count) / real(rate), " seconds"

        !outputファイルへの出力
        write(strN,'(I0)') N                ! 整数
        write(strJ,'(F5.3)') J              ! 小数点以下3桁
        write(strD,'(F5.3)') D
        outputPath = './output/result_' // trim(strN) // '_' // trim(strJ) // '_' // trim(strD) // '_' &
         // vectorPath // '.txt'
        open(unit=10, file=outputPath, status="replace", action="write")
        do L2 = 1,valInd
            write(10, '(I0,",",F12.6,",",I0)') L2,vals(L2),times(L2)
        end do
        close(10)
    end do

    print *, "passed"

contains
    !nCrを計算する
    function nCr(n,r)
        integer, intent(in) :: n,r
        integer :: nCr
        integer :: i
        nCr = 1
        do i = 1, r
            nCr = nCr * (n-i+1) / (i)
        end do
    end function nCr

    !コンビネーション全通りを上げる
    function combinations(n,r) result(combList)
        integer, intent(in) :: n,r
        integer, allocatable :: combList(:)
        integer:: L1
        integer :: counter,createCounter = 1,bit
        integer, allocatable :: tempConb(:)
        allocate(tempConb(n))
        allocate(combList(nCr(n,r)))
        createCounter = 1
        do L1 = 0,2**n-1
            counter = 0
            do L2 = 0,n-1
                bit = getbit(L1,L2)
                counter = counter + bit
                tempConb(L2+1) = bit
                if(counter > r) then
                    exit
                end if
            end do
            if(counter == r)then
                combList(createCounter) = L1
                createCounter = createCounter + 1
            end if
        end do
    end function combinations

    !ある値を二進数にしたときのrank番目を出力
    function getbit(num,rank) result(bit)
        integer, intent(in) :: num
        integer,intent(in)::rank
        integer::bit
        bit = mod((num/(2**rank)),2)
    end function getbit

    !ある値のrank番目にinsertBitを挿入した値を返却
    function setbit(num,rank,insertBit) result(newnum)
        integer, intent(in) :: num,rank,insertBit
        integer::newnum
        newnum = (num/(2**rank)*2+insertBit)*2**rank+mod(num,(2**rank))
    end function setbit

    !ある状態ベクトルを与えられたらH*xを返す
    function Hamil(vecData,Bonds,J,D,ALL_STATE_NUM,UP_NUM,restStates,restStatesNum,idToIndChash) result(retHamil)
        integer, intent(in)::ALL_STATE_NUM,restStatesNum,restStates(restStatesNum),idToIndChash(0:2**N-1)
        complex, intent(in)::vecData(ALL_STATE_NUM)
        integer, intent(in)::Bonds(bond_Num,3),UP_NUM
        real, intent(in)::J,D
        complex::retHamil(ALL_STATE_NUM)
        integer, allocatable :: combList(:)
        integer::L1,L2,L3,beforeSite,afterSite,mod
        integer::state,fripstate,stateInd,fripstateInd
        do L1 = 1,ALL_STATE_NUM
            retHamil(L1) = J/4.0*vecData(L1)*bond_Num
        end do
        do L1 = 1,bond_Num
            !各ボンドについて↑↓と↓↑について調べたい
            beforeSite = Bonds(L1,1)
            afterSite = Bonds(L1,2)
            mod = Bonds(L1,3)
            do L2 = 1,restStatesNum
                if(afterSite < beforeSite) then
                    state = setbit(setbit(restStates(L2),afterSite,0),beforeSite,1)
                    fripstate = setbit(setbit(restStates(L2),afterSite,1),beforeSite,0)
                else
                    state = setbit(setbit(restStates(L2),beforeSite,1),afterSite,0)
                    fripstate = setbit(setbit(restStates(L2),beforeSite,0),afterSite,1)
                end if
                stateInd = idToIndChash(state)
                fripStateInd = idToIndChash(fripstate)
                retHamil(stateInd)     = retHamil(stateInd)-0.5 * J * vecData(stateInd)
                retHamil(fripStateInd) = retHamil(fripStateInd)-0.5 * J * vecData(fripStateInd)
                retHamil(stateInd)    = retHamil(stateInd) + (cmplx(0.5 * J , +0.5 *mod*D)) * vecData(fripStateInd)
                retHamil(fripStateInd) = retHamil(fripStateInd) + (cmplx(0.5 * J ,- 0.5 * mod * D)) * vecData(stateInd)
            end do
        end do
    end function Hamil

    !与えられた状態番号がStates配列の何番目かを返す関数
    function idToInd(StateNum,SiteNum,UP_NUM) result(PattarnCount)
        integer, intent(in)::UP_NUM,SiteNum
        integer::StateNum
        integer::L1,bitcount,bit
        integer::PattarnCount
        bitcount = 0
        PattarnCount = 1
        do L1 = 0,SiteNum-1
            bitcount = bitcount + getbit(StateNum,L1)
        end do
        if(bitcount /= UP_NUM)then
            print *,"invaild state"
            PattarnCount = -1
        else
            bitcount = 0
            do L1 = 1,SiteNum
                bit = getbit(StateNum,SiteNum-L1)
                if(bit == 1)then
                    PattarnCount = PattarnCount + nCr(SiteNum-L1,UP_NUM-bitcount)
                    bitcount = bitcount+1
                end if
            end do
        end if
    end function idToInd

    !ランダムな正規化された複素配列を返す
    function random_vec_gene(leng) result(vector)
        integer, intent(in)::leng
        complex::vector(leng)
        integer::L1
        real::a,b,norm

        call random_seed()
        norm = 0
        do L1 = 1,leng
            call random_number(a)
            call random_number(b)
            vector(L1) = cmplx(a - 0.5d0, b - 0.5d0, kind=8)
        end do
        do L1 = 1,leng
            norm = norm + abs(vector(L1))**2
        end do
        vector = vector / sqrt(norm)
    end function random_vec_gene

    !レイリー商の自作
    function raily(vecData,veclength,accuracy,Bonds,J,D,UP_NUM,restStates,restStatesNum,chashPath,idToIndChash) result(val)
        integer,intent(in)::veclength,Bonds(bond_Num,3),UP_NUM,restStatesNum,restStates(restStatesNum),idToIndChash(0:2**N-1)
        complex,intent(in)::vecData(veclength)
        real,intent(in)::accuracy,J,D

        integer::L1,maxLoops = 10000,L2
        complex::newVec(veclength),beforeVec(veclength)
        real::val,r,attenuation,norm,rand,beforeR

        !学習率計算用の変数
        real::X_HX,dX_HX,X_HdX,dX_HdX,X_X,X_dX,dX_dX
        real::maxAttenution,minAttenution,firstAttenution,secondAttenution
        real::firstR,secondR

        !ファイル書き込み用
        character(len=50) :: strN, strJ, strD, strU,strReal,strIma
        character(len=50),intent(in)::chashPath

        beforeVec(:) = vecData(:)
        !デフォルトの学習率
        attenuation = 0.15/(log10(veclength*1.0)*0.6)
        maxLoops = 3000
        r = 0
        beforeR = r
        val = 100
        Do L1 = 1,maxLoops
            !HXを計算する
            newVec = Hamil(beforeVec,Bonds,J,D,veclength,UP_NUM,restStates,restStatesNum,idToIndChash)
            !収束条件用に一つ前のエネルギーを保管しておく
            beforeR = r
            r = 0
            !レイリー商の計算
            Do L2 = 1,veclength
                r = r + conjg(beforeVec(L2))*newVec(L2)
            end Do
            norm = 0
            !残差ベクトルを作成
            Do L2 = 1,veclength
                norm = norm + abs(newVec(L2)-r*beforeVec(L2))**2
            end Do
            norm = sqrt(norm)
            !収束しているなら抜ける(L1でたまたまr=0近傍の場合を回避/たまにある極地もどきを回避)
            if(abs(r-beforeR) < accuracy .and. L1 /= 1 .and. norm < 0.05)then
                print *,L1
                val = r
                exit
            end if 
            !収束してない場合たまにaccuracyの何倍かを教えるし、学習率も更新する
            call random_number(rand) 
            if(rand < veclength/10.0**7.0)then
                print *,""
                print *,"Current UP_NUM = ",UP_NUM
                print *,"number",L1
                print *,"conver is ",abs(r-beforeR) / accuracy," times"
                print *,"Enegy is",r
                print *,"norm ",norm
            end if 
            !ベクトルを更新し正規化
            norm = 0
            do L2 = 1,veclength
                beforeVec(L2) = beforeVec(L2) - attenuation*(newVec(L2)-r*beforeVec(L2))
                norm = norm + abs(beforeVec(L2))**2
            end do
            !正規化
            do L2 = 1,veclength
                beforeVec(L2) = beforeVec(L2) / sqrt(norm)
            end do
        end do
        !周回回数が巨大な場合は落として、現在の値を返却
        if(L1-1 == maxLoops) then
            print *,"not converge",attenuation
            val = r
        end if
    end function raily
end program baseEnegy
