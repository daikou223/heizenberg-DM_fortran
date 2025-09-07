program baseEnegy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    !可変部
    real,parameter::J =1.,D = 0.1,accuracy = 1e-5
    integer,parameter::N = 24,maxSearch = 1
    integer,parameter::bond_Num = 2*n
    complex,parameter::reaser = 30
    real,parameter::debugNum = 0.01

    !非可変部
    character(len = 100)::inputPath,outputPath,chashPath,statePath,vectorPath
    character(len = 100)::tmpChar,tmpRead1,tmpRead2,tmpRead
    character(len = 10)::vectors(maxSearch,2,2)
    integer,allocatable::states(:)
    integer::Bonds(bond_Num,3)
    integer::ALL_STATE_NUM,valInd
    real::enegys(N),val,vals(N)
    integer:: v1,t1,t2,t3,t4,flag,restStatesNum
    integer::l1,l2,UP_NUM,DOWN_NUM,idToIndChach(0:2**N-1)
    integer, allocatable :: combList(:),restStates(:)
    complex, allocatable::newState(:),initState(:)
    integer::testnum
    complex,allocatable::ret_hamil(:),initData(:),testVec(:),testInitVec(:)
    integer :: start_count, end_count, rate, count_max, times(N),before_count
    !ファイル書き込み用
    character(len=50) :: strN, strJ, strD

    write(tmpChar,'(I2)') N
    inputPath = ".\input\N" // trim(adjustl(tmpchar)) // "\List.txt"
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
    do l1 = 1,maxSearch
        call system_clock(start_count, rate, count_max)
        before_count = start_count
        flag = 1
        vectorPath = "(" // trim(adjustl(vectors(l1,1,1))) // &
                     "," // trim(adjustl(vectors(l1,1,2))) // &
                     ")_(" // trim(adjustl(vectors(l1,2,1))) // &
                     "," // trim(adjustl(vectors(l1,2,2))) // &
                     ")" // "_N" // trim(adjustl(tmpchar)) // ".txt"
        statePath = ".\input\N" // trim(adjustl(tmpchar)) // &
                     "\" // vectorPath 
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
                states = combinations(N,UP_NUM)
                restStates = combinations(N-2,UP_NUM-1)
                restStatesNum = nCr(N-2,UP_NUM-1)
                print *,"create states"
                do L2 = 1,ALL_STATE_NUM
                    idToIndChach(states(L2)) = idToInd(states(L2),N,UP_NUM)
                end do
                print *,"chashed IdToInd"
                initState = random_vec_gene(ALL_STATE_NUM)
                print *,"start raily method"
                val = raily(initState,ALL_STATE_NUM,accuracy,Bonds,J,D,UP_NUM,idToIndChach,restStates,restStatesNum)
                print *,val
                vals(UP_NUM) = val
                valInd = UP_NUM
                call system_clock(end_count)
                times(UP_NUM) = real(end_count - before_count) / real(rate)
                beforeCount = end_count
            end if
        end do

        call system_clock(end_count)
        print *, "Execution time (wall): ", real(end_count - start_count) / real(rate), " seconds"

        !outputファイルへの出力
        write(strN,'(I0)') N                ! 整数
        write(strJ,'(F8.3)') J              ! 小数点以下3桁
        write(strD,'(F8.3)') D
        outputPath = '.\output\result_' // trim(strN) // '_' // trim(strJ) // '_' // trim(strD) // '_' &
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

    !最小固有値を出さないといけないので、適当な大きい値を足して、-をひっかける
    function HamilWarpper(vecData,Bonds,J,D,ALL_STATE_NUM,UP_NUM,idToIndChach,restStates,restStatesNum) result(retHamil)
        integer, intent(in)::ALL_STATE_NUM,idToIndChach(0:2**N-1),restStatesNum,restStates(restStatesNum)
        complex, intent(in)::vecData(ALL_STATE_NUM)
        integer, intent(in)::Bonds(bond_Num,3),UP_NUM
        real, intent(in)::J,D
        complex::retHamil(ALL_STATE_NUM)
        integer::L1
        retHamil = Hamil(vecData,Bonds,J,D,ALL_STATE_NUM,UP_NUM,idToIndChach,restStates,restStatesNum)
        do L1 = 1,ALL_STATE_NUM
            retHamil(L1) = -1*retHamil(L1) + reaser*vecData(L1)
        end do
    end function HamilWarpper
    !ある状態ベクトルを与えられたらH*xを返す
    function Hamil(vecData,Bonds,J,D,ALL_STATE_NUM,UP_NUM,idToIndChach,restStates,restStatesNum) result(retHamil)
        integer, intent(in)::ALL_STATE_NUM,idToIndChach(0:2**N-1),restStatesNum,restStates(restStatesNum)
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
                stateInd = idToIndChach(state)
                fripStateInd = idToIndChach(fripstate)
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

    function beki(vecData,veclength,accuracy,Bonds,J,D,UP_NUM,idToIndChach,restStates,restStatesNum) result(val)
        integer,intent(in)::veclength,Bonds(bond_Num,3),UP_NUM,idToIndChach(0:2**N-1),restStatesNum,restStates(restStatesNum)
        complex,intent(in)::vecData(veclength)
        real,intent(in)::accuracy,J,D

        integer::L1,maxLoops = 10000,L2
        complex::newVec(veclength),normalNewVec(veclength),beforeVec(veclength)
        real::val,difNorm,newVecNorm

        beforeVec(:) = vecData(:) 
        do L2 = 1,maxLoops
            newVec = HamilWarpper(beforeVec,Bonds,J,D,ALL_STATE_NUM,UP_NUM,idToIndChach,restStates,restStatesNum)
            newVecNorm = 0
            do L1 = 1,veclength
                newVecNorm = newVecNorm  + abs(newVec(L1))**2
            end do
            newVecNorm  = sqrt(newVecNorm)
            do L1 = 1,veclength
                normalNewVec(L1) = newVec(L1)/newVecNorm
            end do
            difNorm = 0
            do L1 = 1,veclength
                difNorm = difNorm + abs(beforeVec(L1)-normalNewVec(L1))**2
            end do
            difNorm = sqrt(difNorm)
            call random_number(r)  
            if (difNorm < accuracy) then
                val = ((newVec(1)/beforeVec(1))-reaser)*(-1)
                exit
            else
                beforeVec(:) = normalNewVec(:)
            end if
        end do
        if(L2 == maxLoops) then
            val = -100000.0
            print *,"not converge"
        end if
    end function beki

    !レイリー商の自作
    function raily(vecData,veclength,accuracy,Bonds,J,D,UP_NUM,idToIndChach,restStates,restStatesNum) result(val)
        integer,intent(in)::veclength,Bonds(bond_Num,3),UP_NUM,idToIndChach(0:2**N-1),restStatesNum,restStates(restStatesNum)
        complex,intent(in)::vecData(veclength)
        real,intent(in)::accuracy,J,D

        integer::L1,maxLoops = 10000,L2
        complex::newVec(veclength),normalNewVec(veclength),beforeVec(veclength),nablaR(veclength),HdX(veclength)
        real::val,r,attenuation,norm,rand

        !学習率計算用の変数
        real::X_HX,dX_HX,X_HdX,dX_HdX,X_X,X_dX,dX_dX
        real::maxAttenution,minAttenution,firstAttenution,secondAttenution
        real::firstR,secondR

        beforeVec(:) = vecData(:)
        !デフォルトの学習率
        attenuation = 0.15
        maxLoops = 100000
        r = 0
        val = 100
        Do L1 = 1,maxLoops
            newVec = Hamil(beforeVec,Bonds,J,D,ALL_STATE_NUM,UP_NUM,idToIndChach,restStates,restStatesNum)
            r = 0
            Do L2 = 1,veclength
                r = r + conjg(beforeVec(L2))*newVec(L2)
            end Do
            norm = 0
            Do L2 = 1,veclength
                nablaR(L2) = (newVec(L2)-r*beforeVec(L2))
                norm = norm + abs(nablaR(L2))**2
            end Do
            call random_number(rand) 
            if(norm < accuracy)then
                print *,L1
                val = r
                exit
            end if 
            if(rand < 0) then
                attenuation = -1
            end if 
            if(attenuation < 0) then
                HdX = HamilWarpper(nablaR,Bonds,J,D,ALL_STATE_NUM,UP_NUM,idToIndChach,restStates,restStatesNum) 
                X_HX = r
                dX_HX = 0
                X_HdX = 0
                dX_HdX = 0
                X_X = 1
                X_dX = 0
                dX_X = 0
                dX_dX = 0
                do L2 = 1,veclength
                    X_HX = X_HX + conjg(beforeVec(L2))*newVec(L2)
                    dX_HX = dX_HX + conjg(nablaR(L2))*newVec(L2)
                    X_HdX = X_HdX + conjg(beforeVec(L2))*HdX(L2)
                    dX_HdX = dX_HdX + conjg(nablaR(L2))*HdX(L2)
                    X_dX = X_dX + conjg(beforeVec(L2))*nablaR(L2)
                    dX_X = dX_X + conjg(nablaR(L2))*beforeVec(L2)
                    dX_dX = dX_dX + conjg(nablaR(L2))*nablaR(L2)
                end do
                maxAttenution = 1
                minAttenution = 0
                do L2 = 1,40
                    firstAttenution = (2*minAttenution + maxAttenution)/3.0
                    secondAttenution = (minAttenution + 2*maxAttenution)/3.0
                    firstR = (X_HX-firstAttenution*dX_HX-firstAttenution*X_HdX+(firstAttenution**2)*dX_HdX)/&
                        (X_X-firstAttenution*X_dX-firstAttenution*dX_X+(firstAttenution**2)*dX_dX)
                    secondR = (X_HX-secondAttenution*dX_HX-secondAttenution*X_HdX+(secondAttenution**2)*dX_HdX)/&
                        (X_X-secondAttenution*X_dX-secondAttenution*dX_X+(secondAttenution**2)*dX_dX)
                    if(firstR < secondR) then
                        maxAttenution = secondAttenution
                    else
                        minAttenution = firstAttenution
                    end if
                end do
                attenuation = maxAttenution
            end if
            norm = 0
            do L2 = 1,veclength
                beforeVec(L2) = beforeVec(L2) - attenuation*nablaR(L2)
                norm = norm + abs(beforeVec(L2))**2
            end do
            do L2 = 1,veclength
                beforeVec(L2) = beforeVec(L2) / sqrt(norm)
            end do
        end do
        if(L1-1 == maxLoops) then
            print *,"not converge",attenuation
            val = r
        end if
    end function raily
end program baseEnegy