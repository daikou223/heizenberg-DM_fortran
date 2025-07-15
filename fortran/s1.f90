!************* Sample program #1 ***************
!   Eigenvalues and an eigenvector / lnc1, lncv1
!    Precision check and correlation functions
!***********************************************
      parameter (n=21,idim=352716,ibond=n*2)
      implicit real*8 (a-h,o-z)
      dimension E(4)
      dimension list1(idim),list2(2,0:2**15)
      dimension bondwt(ibond),ipair(2*ibond),zrtio(ibond)
      dimension npair(2)
      dimension wk(idim,2)
      dimension x(idim)

      data bondwt/ibond*-1.0d0/
      data zrtio/ibond*1.0d0/
      data ipair/1,21, 1,8,&
      & 2,6, 2,21, &
      & 3,2, 3,20,&
      & 4,3, 4,7,&
      & 5,4, 5,1,&
      & 6,3, 6,10,&
      & 7,5, 7,12,&
      & 8,5, 8,13,&
      & 9,6, 9,8,&
      & 10,9, 10,14,&
      & 11,10, 11,7,&
      & 12,11, 12,19,&
      & 13,9, 13,16,&
      & 14,11, 14,18,&
      & 15,13, 15,12,&
      & 16,15, 16,20,&
      & 17,16, 17,14,&
      & 18,17, 18,1,&
      & 19,15, 19,2,&
      & 20,17, 20,4,&
      & 21,19, 21,18/
      nvec=1
      iv=idim/3
      call sz(n,idim,0.5d0,list1,list2)

!*** Eigenvalues
      call lnc1(n,idim,ipair,bondwt,zrtio,ibond,&
     &                 nvec,iv,E,itr,wk,idim,list1,list2)
      print 100,e,itr
 100  format(/' [Eigenvalues]  '/2x,4f14.8&
     &       /' [Iteration number]'/i8)

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
      end

