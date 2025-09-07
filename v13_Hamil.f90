! Hamil.f90
module hamil_mod
  implicit none
contains

  subroutine Hamil_fortran(vecData, ret_hamil, Bonds, numBonds, UP_NUM, allComb, nStates, idToInd)
    implicit none
    integer, intent(in) :: numBonds, UP_NUM, nStates
    integer, intent(in) :: idToInd(0:nStates-1)
    integer, intent(in) :: Bonds(numBonds,3)
    integer, intent(in) :: allComb(:,:,:)
    complex(8), intent(in) :: vecData(nStates)
    complex(8), intent(out) :: ret_hamil(nStates)

    integer :: bondInd, combInd, stateIndex, fripStateIndex
    integer :: i, upSite, combSize
    integer :: state, fripState

    ! 初期化
    ret_hamil = (0.0d0,0.0d0)

    ! 定数 J
    real(8) :: bond_mod_J
    bond_mod_J = 1.0d0

    ! ベースの対角成分
    do i = 1, numBonds
       ret_hamil = ret_hamil + (0.25d0 * bond_mod_J * numBonds) * vecData
    end do

    ! 各ボンドの操作
    do bondInd = 1, numBonds
       combSize = size(allComb(bondInd,:,:),2)
       do combInd = 1, combSize
          ! 状態生成
          state = 2**(Bonds(bondInd,1)-1)
          fripState = 2**(Bonds(bondInd,2)-1)
          do i = 1, UP_NUM-1
             upSite = allComb(bondInd,combInd,i)
             state = ior(state, 2**(upSite-1))
             fripState = ior(fripState, 2**(upSite-1))
          end do
          stateIndex = idToInd(state)
          fripStateIndex = idToInd(fripState)

          ! ハミルトニアン作用
          ret_hamil(stateIndex)     = ret_hamil(stateIndex)     - 0.5d0 * bond_mod_J * vecData(stateIndex)
          ret_hamil(fripStateIndex) = ret_hamil(fripStateIndex) - 0.5d0 * bond_mod_J * vecData(fripStateIndex)
          ret_hamil(fripStateIndex) = ret_hamil(fripStateIndex) + (0.5d0*bond_mod_J - 0.5d0*cmplx(0.0d0,Bonds(bondInd,3))) * vecData(stateIndex)
          ret_hamil(stateIndex)     = ret_hamil(stateIndex)     + (0.5d0*bond_mod_J + 0.5d0*cmplx(0.0d0,Bonds(bondInd,3))) * vecData(fripStateIndex)
       end do
    end do

  end subroutine Hamil_fortran

end module hamil_mod
