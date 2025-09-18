program baseEnegy
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none


    !可変部
    real,parameter::J =1.,D = 0.1,accuracy = 1e-5
    integer,parameter::N = 6,maxSearch = 1
    integer,parameter::bond_Num = 2*N

    !非可変部-file
    character(len = 100)::inputPath,outputPath,chashPath,statePath,vectorPath
    character(len = 100)::tmpChar,tmpRead1,tmpRead2,tmpRead,ReadChash
    character(len = 10)::vectors(maxSearch,2,2)
    integer:: t1,t2,t3,t4,beforeSite,afterSite,mod
    character(len=50) :: strN, strJ, strD, strU
    logical::fileExists

    !非可変部-main
    integer,allocatable::states(:)
    integer::Bonds(bond_Num,3)
    integer::ALL_STATE_NUM,valInd
    real::val,vals(N),vecDataNorm
    integer:: restStatesNum
    integer, allocatable :: combList(:),restStates(:)
    complex, allocatable::newState(:),vecData(:)
    complex,allocatable::ret_hamil(:),initData(:)
    real::railyDivs(N)

    !非可変部-raily
    real::attenutation
    integer,parameter::maxloops = 3000
    integer::railyLoop,restStatesInd
    complex,allocatable::newVecData(:)
    integer::maxSite,minSite
    integer::minStateInd,maxStateInd,predictStateInd
    integer::stateInd,fripStateInd
    real::railyDiv,nablaRNorm,nextVecDataNorm,attenuation,innerRate
    integer:: state, fripState, maxMask, minMask, restMask
    integer(8):: maxMask_, minMask_, restMask_
    integer::mask,key,upperUseKey

    !非可変部-乱数
    real :: a,b,rand

    !非可変部-chash
    integer::nCr(0:N,0:N)
    integer,allocatable::idToInd(:)
    integer,parameter::collapseList(55) =  [ 527744199 , 956575440, 496739444,670680204,111485899,530512614,&
     147420730, 970446391, 331698202, 571882804, 242025422, 864200478, 961745359, 755734019, &
     503753589, 381178730, 902675107, 351434920, 489318112, 539076824, 294007524, 34155980, &
     692702580, 601487816, 623492890, 686261039, 883346590, 1066239057, 1054396908, 621034762,&
      1037879269, 196138108, 979317313, 133899874, 206149928, 954313712, 74069865, 90810165, &
      700190095, 109179671, 843423738, 526934554, 425697326, 444232068, 545474910, 684040096,&
       1027339630, 1006583862, 690409966, 631407572, 537654983, 158950028, 566370613, 628623012, 412664663]

    !非可変部-tmp変数
    integer::tmpnCr,tmpBeforeStatesCounter,tmpStatesCounter
    integer,allocatable::tmpStates(:,:),tmpBeforeStates(:,:)

    !非可変部-loop
    integer::l1,l2,UP_NUM,DOWN_NUM,vecInd,bondInd,tmpBeforeStatesInd

    !非可変部-time
    integer :: start_count, end_count, rate, count_max, times(N)

    !非可変部-debug
    integer::testnum,beforeInd,afterInd
    complex,allocatable::testVec(:),testInitVec(:)

    !処理はじめ**************************************
    !周期境界条件List.txtから取得
    allocate(idToInd(0:2**(N-2)-1))
    write(tmpChar,'(I2)') N
    inputPath = "./input/N" // trim(adjustl(tmpchar)) // "/List.txt"
    open(unit=10, file=inputPath, status="old")
        do vecInd = 1,maxSearch
            read(10,'(A)') tmpRead1
            read(tmpRead1,*) t1,t2,t3,t4
            write(vectors(vecInd,1,1),'(I5)') t1
            write(vectors(vecInd,1,2),'(I5)') t2
            write(vectors(vecInd,2,1),'(I5)') t3
            write(vectors(vecInd,2,2),'(I5)') t4
        end do
    close(10)
    !上記で取得できたものうちいい奴に対して計算を行う
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
            do bondInd = 1,bond_Num
                read(10,'(A)') tmpRead1
                read(tmpRead1,*) beforeSite,afterSite,mod
                Bonds(bondInd,1) = beforeSite
                Bonds(bondInd,2) = afterSite
                Bonds(bondInd,3) = mod
            end do
        close(10)
        !nCrを作る
        nCr(0,0) = 1
        do l1 = 1,N
            nCr(l1,0) = 1
            nCr(l1,l1) = 1
            do l2 = 1,l1-1
                nCr(l1,l2) = nCr(l1-1,l2-1) + nCr(l1-1,l2)
            end do
        end do
        do UP_NUM = N,1,-1
            print *,"UP_NUM = ",UP_NUM
            DOWN_NUM = N-UP_NUM
            !逆のときは対称性より計算しない
            if(UP_NUM <= DOWN_NUM) then
                !ベクトルの大きさを作る
                ALL_STATE_NUM = nCr(N,UP_NUM)
                print *,"HamilSize ",ALL_STATE_NUM
                !全状態を作る
                !class的なもの(state,一の数)を持つ
                if (allocated(tmpStates)) deallocate(tmpStates)
                allocate(tmpStates(ALL_STATE_NUM,2))
                if (allocated(tmpBeforeStates)) deallocate(tmpBeforeStates)
                allocate(tmpBeforeStates(ALL_STATE_NUM,2))
                if (allocated(states)) deallocate(states)
                allocate(states(ALL_STATE_NUM))
                !幅優先を再起じゃなく書く
                tmpBeforeStates(1,:) = [0, 0]
                tmpBeforeStatesCounter = 1
                do l1 = 1,N
                    !それぞれに対して、上限下限に引っかかってないなら
                    !1,0を追加したもの追加
                    tmpStatesCounter = 0 !tmpStateのサイズを持っておく
                    do tmpBeforeStatesInd = 1,tmpBeforeStatesCounter
                        if(tmpBeforeStates(tmpBeforeStatesInd,2) > l1 +UP_NUM -N -1) then 
                            !末尾に0を足す(つまりは2倍)で1の数は変わらない
                            tmpStatesCounter = tmpStatesCounter + 1
                            tmpStates(tmpStatesCounter,1) = tmpBeforeStates(tmpBeforeStatesInd,1)*2
                            tmpStates(tmpStatesCounter,2) = tmpBeforeStates(tmpBeforeStatesInd,2)
                        end if
                        if(tmpBeforeStates(tmpBeforeStatesInd,2)  < UP_NUM) then 
                            !末尾に1を足す
                            tmpStatesCounter = tmpStatesCounter + 1
                            tmpStates(tmpStatesCounter,1) = tmpBeforeStates(tmpBeforeStatesInd,1)*2+1
                            tmpStates(tmpStatesCounter,2) = tmpBeforeStates(tmpBeforeStatesInd,2)+1
                        end if
                    end do
                    tmpBeforeStates = tmpStates
                    tmpBeforeStatesCounter = tmpStatesCounter
                end do
                idToInd(:) = -1
                upperUseKey = 0
                do l1 = 1,ALL_STATE_NUM
                    states(l1) = tmpStates(l1,1) 
                    key = ishft(states(l1),-2)
                    mask = ishft(1,N-2)-1
                    do L2 = 1,size(collapseList)
                        if(idToInd(key) > 0) then
                            key = ieor(ishft(states(l1),-2),iand(mask,collapseList(L2)))
                            if(l2 > upperUseKey) then
                                print *,"DangerousUseKey:",l2
                                upperUseKey = l2
                            end if
                        else 
                            exit
                        end if
                    end do
                    if (l2-1 == size(collapseList))then
                        print *,"key wither"
                        print *,testvec(100)
                    end if
                    idToInd(key) = l1
                end do
                !テスト*********************************************************************
                ! print*,"test"
                ! mask = (ishft(1,N))-1
                ! do l1 = 1,ALL_STATE_NUM
                !     print *,states(l1)
                !     if(idToInd(ishft(states(l1),-2)) < ALL_STATE_NUM) then
                !         print *,"stetas(",idtoind(ishft(states(l1),-2)),") = ",states(idtoind(ishft(states(l1),-2))),"1"
                !     else
                !         beforeInd =  ishft(idToInd(ishft(states(l1),-2)),-N)
                !         afterInd =  iand(idToInd(ishft(states(l1),-2)),mask)
                !         if(states(beforeInd) == states(l1)) then
                !             print *,"stetas(",beforeInd,") = ",states(beforeInd),"2"
                !         else
                !             print *,"stetas(",afterInd,") = ",states(afterInd),"3"
                !         end if
                !     end if
                ! end do
                !********************************************************
                print *,"created state"
                !残状態を作る
                restStatesNum = nCr(N-2,UP_NUM-1)
                !class的なもの(state,一の数)を持つ
                if (allocated(tmpStates)) deallocate(tmpStates)
                allocate(tmpStates(restStatesNum,2))
                if (allocated(tmpBeforeStates)) deallocate(tmpBeforeStates)
                allocate(tmpBeforeStates(restStatesNum,2))
                if (allocated(restStates)) deallocate(restStates)
                allocate(restStates(restStatesNum))
                !幅優先を再起じゃなく書く
                tmpBeforeStates(1,:) = [0, 0]
                tmpBeforeStatesCounter = 1
                do l1 = 1,N-2
                    !それぞれに対して、上限下限に引っかかってないなら
                    !1,0を追加したもの追加
                    tmpStatesCounter = 0 !tmpStateのサイズを持っておく
                    do tmpBeforeStatesInd = 1,tmpBeforeStatesCounter
                        if(tmpStatesCounter+1 > restStatesNum)then
                            print *,"ここ"
                        end if
                        if(tmpBeforeStates(tmpBeforeStatesInd,2) > l1 +(UP_NUM-1) -(N-2) -1) then 
                            !末尾に0を足す(つまりは2倍)で1の数は変わらない
                            tmpStatesCounter = tmpStatesCounter + 1
                            tmpStates(tmpStatesCounter,1) = tmpBeforeStates(tmpBeforeStatesInd,1)*2
                            tmpStates(tmpStatesCounter,2) = tmpBeforeStates(tmpBeforeStatesInd,2)
                        end if
                        if(tmpBeforeStates(tmpBeforeStatesInd,2)  < (UP_NUM-1)) then 
                            !末尾に1を足す
                            tmpStatesCounter = tmpStatesCounter + 1
                            tmpStates(tmpStatesCounter,1) = tmpBeforeStates(tmpBeforeStatesInd,1)*2+1
                            tmpStates(tmpStatesCounter,2) = tmpBeforeStates(tmpBeforeStatesInd,2)+1
                        end if
                    end do
                    tmpBeforeStates = tmpStates
                    tmpBeforeStatesCounter = tmpStatesCounter
                end do

                do l1 = 1,restStatesNum
                    restStates(l1) = tmpStates(l1,1) 
                end do
                print *,"created reststate"
                !ランダムなベクトルを生成
                vecDataNorm = 0
                call random_seed()
                if (allocated(vecData)) deallocate(vecData)
                allocate(vecData(ALL_STATE_NUM))
                do L1 = 1,ALL_STATE_NUM
                    call random_number(a)
                    call random_number(b)
                    vecData(L1) = cmplx(a - 0.5d0, b - 0.5d0, kind=8)
                    vecDataNorm = vecDataNorm + (a-0.5)**2+(b-0.5)**2
                end do
                vecData = vecData/sqrt(vecdataNorm)
                print *,"start raily"
                !レイリー商反復
                attenuation = 0.15/(log10(ALL_STATE_NUM*1.0)*0.6)
                do railyLoop = 1,maxLoops
                    newVecData = (J/4.0) * bond_Num * vecData
                    do bondInd = 1,bond_Num
                        !結合の情報を取得
                        beforeSite = Bonds(bondInd,1)
                        afterSite = Bonds(bondInd,2)
                        mod = Bonds(bondInd,3)
                        !minSiteの影響で一つずれるから
                        maxSite = max(beforeSite,Aftersite)-1
                        minSite = min(beforeSite,Aftersite)
                        !入れ込むためのマスク
                        restMask = (ishft(1,minSite))-1
                        minMask = (ishft(1,maxSite))-1-restMask
                        maxMask = (ishft(1,N))-1-minMask-restmask
                        do restStatesInd = 1,restStatesNum
                            !実際の入れ込み
                            state = iand(restStates(reststatesInd),maxMask)*4 + iand(restStates(reststatesInd),minMask)*2 +&
                            iand(restStates(reststatesInd),restMask)
                            state = ibset(state,beforeSite)
                            fripState = iand(restStates(reststatesInd),maxMask)*4 + iand(restStates(reststatesInd),minMask)*2 +&
                            iand(restStates(reststatesInd),restMask)
                            fripstate = ibset(fripState,afterSite)
                            !まずはstateのほうから
                            mask = (ishft(1,N-2))-1
                            key = (ishft(state,-2))
                            do l1 = 1,size(collapseList)
                                if(States(idToInd(key)) == state) then
                                    exit
                                else
                                    key = ieor(ishft(state,-2),iand(mask,collapseList(L1)))
                                end if
                            end do 
                            stateInd = idToInd(key)
                            !fripstateのほう
                            key = (ishft(fripState,-2))
                            do l1 = 1,size(collapseList)
                                if(States(idToInd(key)) == fripState) then
                                    exit
                                else
                                    key = ieor(ishft(fripState,-2),iand(mask,collapseList(L1)))
                                end if
                            end do 
                            fripStateInd = idToInd(key)
                            if(stateInd < 1 .or. stateInd > ALL_STATE_NUM .or. fripStateInd < 1 .or. &
                             fripStateInd > ALL_STATE_NUM) then
                                print *,stateInd,fripStateInd,state,fripState
                            end if
                            !行列演算をする
                            newVecData(stateInd)     = newVecData(stateInd)-0.5 * J * vecData(stateInd)
                            newVecData(fripStateInd) = newVecData(fripStateInd)-0.5 * J * vecData(fripStateInd)
                            newVecData(stateInd)    = newVecData(stateInd) + &
                                (cmplx(0.5 * J , +0.5 *mod*D)) * vecData(fripStateInd)
                            newVecData(fripStateInd) = newVecData(fripStateInd) + &
                                (cmplx(0.5 * J ,- 0.5 * mod * D)) * vecData(stateInd)
                        end do
                        !ボンド固定、のこりのstateで回し終わり
                    end do
                    !ボンド動かし終わり
                    !レイリー商の計算
                    railyDiv = 0
                    Do l1 = 1,ALL_STATE_NUM
                        railyDiv = railyDiv + conjg(vecData(l1))*newVecData(l1)
                    end Do
                    !残差ベクトルのノルム
                    nablaRNorm = 0
                    Do l1 = 1,ALL_STATE_NUM
                        nablaRNorm = nablaRNorm + abs(newVecData(l1)-railyDiv*vecData(l1))**2
                    end Do
                    nablaRNorm = sqrt(nablaRNorm)
                    !収束してない場合たまにaccuracyの何倍かを教える
                    call random_number(rand) 
                    if(rand < ALL_STATE_NUM/(10.0**7.0) .or. railyLoop == 1)then
                        print *,""
                        print *,"attenution",attenuation
                        print *,"Current UP_NUM = ",UP_NUM
                        print *,"number",railyLoop
                        print *,"Enegy is",railyDiv
                        print *,"norm ",nablaRNorm
                    end if 
                    !収束してたら返す
                    if(nablaRNorm < 0.003)then
                        print*,(railyDiv)
                        print*,(railyLoop)
                        railyDivs(UP_NUM) = railyDiv
                        exit
                    end if
                    !ベクトルを更新し正規化
                    nextVecDataNorm = 0
                    do l1 = 1,ALL_STATE_NUM
                        vecData(l1) = vecData(l1) - attenuation*(newVecData(l1)-railyDiv*vecData(l1))
                        nextVecDataNorm = nextVecDataNorm + abs(vecData(l1))**2
                    end do
                    !正規化
                    vecData = vecData / sqrt(nextVecDataNorm)
                end do
                !レイリー商反復の終了
                !最後までいった場合はその時のノルムを返却
                if(L1-1 == maxLoops) then
                    print *,railyLoop
                    railyDivs(UP_NUM) = railyDiv
                    print *,"not converge",nablaRNorm
                end if
            end if
            !UP_NUM <= DOWN_NUMのif
        end do
        print *,"enegys:",railyDivs

        call system_clock(end_count)
        print *, "Execution time (wall): ", real(end_count - start_count) / real(rate), " seconds"

        !outputファイルへの出力
        write(strN,'(I0)') N                ! 整数
        write(strJ,'(F5.3)') J              ! 小数点以下3桁
        write(strD,'(F5.3)') D
        outputPath = './output/result_' // trim(strN) // '_' // trim(strJ) // '_' // trim(strD) // '_' &
         // vectorPath // '.txt'
        open(unit=10, file=outputPath, status="replace", action="write")
        do L2 = 1,N
            write(10, '(I0,",",F12.6,",",I0)') L2,railyDivs(L2)
        end do
        close(10)
    end do
    print *,"passed"
end program baseEnegy
