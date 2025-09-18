program Hamirutonian
    implicit none
    integer,parameter :: particlesNumber = 2
    real(8),parameter :: JAction = 1
    integer situation,bit,newSituation,i,j
    real(8) :: coefficient
    real(8),dimension(2**particlesNumber) :: particlesBits,transition = 0,spins
    real(8),dimension(0:2**particlesNumber-1,0:2**particlesNumber-1) :: reprensentationMatrix = 0.0d0
    integer,dimension(particlesNumber) :: situationBits,copySituationBits
    !全てのスピン状態(2**Nパターン)について調べる
    print *,2**particlesNumber-1
    i = 0
    do situation = 0,2**particlesNumber-1
        !スピン状態を2進数にする
        situationBits = NumberToBit(situation)
        !周期境界以外の各ビットについて行う
        do bit = 1,particlesNumber-1
            !スピン反転するのでコピーしておく
            copySituationBits = situationBits
            !スピン方向が同じならS+orS-の影響で0になる
            !スピン方向が|01>なら
            if(situationBits(bit) == 0 .AND. situationBits(bit+1) == 1)then
                copySituationBits(bit) = 1
                copySituationBits(bit+1) = 0
                newSituation = BitToNumber(copySituationBits)
                print *,29,newSituation
                reprensentationMatrix(situation,newSituation) = reprensentationMatrix(situation,newSituation) + JAction/2
            !スピン方向が|10>なら
            else if(situationBits(bit) == 1 .AND. situationBits(bit+1) == 0)then
                copySituationBits(bit) = 0
                copySituationBits(bit+1) = 1
                newSituation = BitToNumber(copySituationBits)
                print *,36,newSituation
                reprensentationMatrix(situation,newSituation) = reprensentationMatrix(situation,newSituation) + JAction/2
            end if
        end do
        if(particlesNumber /= 2)then
            !境界条件部分
            copySituationBits = situationBits
            !スピン方向が同じならS+orS-の影響で0になる
            !スピン方向が|0...1>なら
            if(situationBits(1) == 0 .AND. situationBits(particlesNumber) == 1)then
                copySituationBits(particlesNumber) = 0
                copySituationBits(1) = 1
                newSituation = BitToNumber(copySituationBits)
                print *,48,newSituation
                reprensentationMatrix(situation,newSituation) = reprensentationMatrix(situation,newSituation) + JAction/2.0
            !スピン方向が|1...0>なら
            else if(situationBits(particlesNumber) == 1 .AND. situationBits(1) == 0)then
                copySituationBits(particlesNumber) = 0
                copySituationBits(1) = 1
                newSituation = BitToNumber(copySituationBits)
                print *,55,newSituation
                reprensentationMatrix(situation,newSituation) = reprensentationMatrix(situation,newSituation) + JAction/2.0
            end if
        end if
        !対角成分
        do bit = 1,particlesNumber-1
            !スピン反転するのでコピーしておく
            copySituationBits = situationBits
            !スピン方向が異なるなら
            if(situationBits(bit) /= situationBits(bit+1))then
                reprensentationMatrix(situation,situation) = reprensentationMatrix(situation,situation) - JAction/4.0
            !スピン方向が同じなら
            else if(situationBits(bit) == situationBits(bit+1))then
                reprensentationMatrix(situation,situation) = reprensentationMatrix(situation,situation) + JAction/4.0
            end if
        end do
        if(particlesNumber /= 2)then
            !境界条件
            !スピン方向が異なるなら
            if(situationBits(1) /= situationBits(particlesNumber))then
                reprensentationMatrix(situation,situation) = reprensentationMatrix(situation,situation) - JAction/4.0
            !スピン方向が同じなら
            else if(situationBits(particlesNumber) == situationBits(1))then
                reprensentationMatrix(situation,situation) = reprensentationMatrix(situation,situation) + JAction/4.0
            end if
        end if
    end do
    do i = 0,2**particlesNumber-1
        print *,"["
        do j = 0,2**particlesNumber-1
            print *,reprensentationMatrix(i,j)
        end do
        print *,"]"
    end do
contains
    function NumberToBit(num)
        integer bit,num,tmp
        integer,dimension(particlesNumber) :: NumberToBit
        tmp = num
        do bit = particlesNumber,1,-1
            NumberToBit(bit) = mod(tmp,2)
            tmp = tmp / 2
        end do
    end function NumberToBit

    function BitToNumber(numBit)
        integer bit,BitToNumber
        integer,dimension(particlesNumber) :: numBit
        BitToNumber = 0
        do bit = 1,size(numBit)
            BitToNumber = BitToNumber + numBit(bit)*2**(size(numBit)-bit)
        end do
    end function BitToNumber
end program Hamirutonian