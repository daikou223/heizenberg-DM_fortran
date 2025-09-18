!************* Sample program #1 ***************
!   Eigenvalues and an eigenvector / lnc1, lncv1
!    Precision check and correlation functions
!***********************************************
program main

parameter (n=28,idim=40116600,ibond=2*n)
implicit real*8 (a-h,o-z)
dimension E(4)
dimension bondwt(ibond),ipair(2*ibond),zrtio(ibond)
dimension npair(2)
integer,allocatable::list1(:),list2(:,:)
real*8,allocatable::wk(:,:),x(:)
real*8::szval=0.d0
integer,allocatable::num_one(:),ncombi(:,:),num_basis(:),num_begin(:)
real::t,secnds
integer::i

allocate(list1(idim),list2(2,0:2**24),wk(idim,2),x(idim))
allocate(num_one(0:2**24-1),ncombi(0:24*2,0:24*2),num_basis(0:2**24-1),num_begin(0:2**24-1))

bondwt(1:n) = -1.0d0
bondwt(n+1:2*n) = -0.5d0
data zrtio/ibond*1.0d0/

!独自に書いたとこ
do i = 1,n-1
      ipair(2*i-1:2*i) = [i,i+1]
end do
ipair(2*n-1:2*n) = [n,1]
do i = 1,n-2
      ipair(2*i-1+2*n : 2*i+2*n) = [i,i+2]
end do
ipair(2*(n-1)-1+2*n: 2*(n-1)+2*n) = [n-1,1]
ipair(2*n-1+2*n : 2*n+2*n) = [n,2]
print *,ipair
!独自終わり

t=secnds(0.0)

      nvec=1
      iv=idim/3
write(*,*) 'sz'
      call szhn(n,idim,szval,list1,list2,num_one,ncombi,num_basis,num_begin)
!      call sz(n,idim,szval,list1,list2)

!*** Eigenvalues
      write(*,*) 'lnc1'
      call lnc1(n,idim,ipair,bondwt,zrtio,ibond,&
     &                 nvec,iv,E,itr,wk,idim,list1,list2)
      print 100,e,itr
 100  format(/' [Eigenvalues]  '/2x,4f14.8&
     &       /' [Iteration number]'/i8)

!stop
!*** Ground-state eigenvector
      call lncv1(n,idim,ipair,bondwt,zrtio,ibond,&
     &            nvec,iv,x,itr,wk,idim,list1,list2)
!- You may alternatively use inv1 / Note: dimension wk(idim,4) -
!      call inv1(n,idim,ipair,bondwt,zrtio,ibond,&
!     &          E(1),iv,x,wk,idim,list1,list2)
!---------------------------------------------------------------
      print *,'[Eigenvector components (selected)]'
      print 120,(x(j),j=13,idim,idim/20)
 120  format(4d18.9)

write(*,*) x(13)

!*** Precision check and correlation functions
      call check1(n,idim,ipair,bondwt,zrtio,ibond,&
     &            x,wk,Hexpec,list1,list2)
      npair(1)=1
      npair(2)=2
      call xcorr(n,idim,npair,1,x,sxx,list1,list2)
      call zcorr(n,idim,npair,1,x,szz,list1)
      print 130,sxx,szz
 130  format(/' [Nearest neighbor correlation functions]'/&
     &       '    sxx :',d18.10,',    szz :',d18.10)

write(*,'("Process Time(sec.):",F10.2)') secnds(t)

end program main

