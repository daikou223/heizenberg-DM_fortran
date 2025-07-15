c***********************************************************
c*                                                         *
c*                    TITPACK Ver. 2                       *
c*                                                         *
c*                    February, 1991                       *
c*                                                         *
c*          Copyright (C) Hidetoshi Nishimori              *
c*                                                         *
c***********************************************************
c
c================== COMMON SUBROUTINES ====================
c
c*** variables marked @ should be given from the main program
c*** variables marked # are evaluated and returned
c
c******** configurations with the specified sz ***********
c
      subroutine szhn(n,idim,szval,list1,list2,num_one,ncombi,
     &                num_basis,num_begin)
c
c    n          @  lattice size
c    idim       @  dimension of the matrix
c    szval      @  total sz
c    list1(i)   #  i-th spin configuration
c    list2      #  inverse list of list1 expressed by the
c                  2-dim search method of M. Ogata and H.Q. Lin.
c
      real*8 szval
      dimension list1(idim),list2(2,0:2**24)
c
      dimension num_one(0:2**24-1)
      dimension ncombi(0:24*2,0:24*2)
      dimension num_basis(0:2**24-1),num_begin(0:2**24-1)
c
      if(szval.lt.-1.0d-13.or.szval.gt.n/2.d0-1.d0+1.d-13)then
        print *,' #(E01)# Variable szval given to sz out of range'
        stop
      end if
      if(idim.lt.3)then
        print *,' #(E02)# Incorrect idim or n given to sz'
        stop
      end if
c
c* initialization
      ihfbit=2**((n+1)/2)
      iupspn=n/2+mod(n,2)+int(szval+0.001d0)
c
c counting one
      do 5000 i=0,2**24-1
        num_one(i)=0
        num_basis(i)=0
        num_begin(i)=0
 5000    continue
      do 5100 i=0,2**24-1
        isz=0
        do 5200 j=0,24-1
          isz=isz+mod(i/2**j,2)
 5200     continue
        num_one(i)=isz 
 5100    continue
c
c number of combination
      do 6000 j=0,24*2
      do 6100 i=0,24*2
        ncombi(i,j)=0
 6100    continue
 6000    continue
      do 6200 i=0,n
        ncombi(i,0)=1
 6200    continue
      do 6500 i=1,n
        do 6400 j=1,i
           if ((i.eq.1).and.(j.eq.1)) then
             ncombi(i,j)=ncombi(i-1,j-1)
           else
             ncombi(i,j)=ncombi(i-1,j-1)+ncombi(i-1,j)
           endif
 6400     continue
 6500     continue
      idim_calc=0
      ib=0
      num_one_ib=num_one(ib)
      if (iupspn.ge.num_one_ib) then
        num_basis(ib)=ncombi((n+1)/2,iupspn-num_one_ib)
      else
        num_basis(ib)=0
      endif
      num_begin(ib)=0
      do 6600 ib=1,2**(n/2)-1
        num_one_ib=num_one(ib)
        if (iupspn.ge.num_one_ib) then
          num_basis(ib)=ncombi((n+1)/2,iupspn-num_one_ib)
        else
          num_basis(ib)=0
        endif
        num_begin(ib)=num_basis(ib-1)+num_begin(ib-1)
 6600    continue
      idim_calc=num_begin(2**(n/2)-1)+num_basis(2**(n/2)-1)
      if (idim_calc.ne.idim) then
        write(*,*) 'Incorrect idim, n or Sz'
      endif
c
c* main loop
!$OMP parallel do 
!$OMP& private(ib,icnt,ibpatn,ja,jb,ia,i,isz),
!$OMP& shared(n,num_begin,ihfbit,num_one,iupspn,idim,list1,list2)
      do 1000 ib=0,2**(n/2)-1
        icnt=num_begin(ib)
        ibpatn=0
        ja=0
        jb=icnt
      do 100 ia=0,2**((n+1)/2)-1
        i=ib*ihfbit+ia    
        isz=num_one(ia)+num_one(ib)
        if(isz.ne.iupspn)go to 100
        icnt=icnt+1
        if(icnt.gt.idim)then
          print *,' #(E02)# Incorrect idim or n given to sz'
          stop
        end if
        if(ib.eq.ibpatn)then
          ja=ja+1
        else
          ibpatn=ib
          ja=1
          jb=icnt-1
        end if
        list1(icnt)=i
        list2(1,ia)=ja
        list2(2,ib)=jb
 100  continue
1000  continue
!$OMP end parallel do 
c
c      if(icnt.eq.idim)return
c      print *,' #(E02)# Incorrect idim or n given to sz'
c      stop
      return
      end
c
c******** configurations with the specified sz ***********
c
      subroutine sz(n,idim,szval,list1,list2)
c
c    n          @  lattice size
c    idim       @  dimension of the matrix
c    szval      @  total sz
c    list1(i)   #  i-th spin configuration
c    list2      #  inverse list of list1 expressed by the
c                  2-dim search method of M. Ogata and H.Q. Lin.
c
      real*8 szval
c      integer*8 i
      dimension list1(idim),list2(2,0:2**24)
c
      if(szval.lt.-1.0d-13.or.szval.gt.n/2.d0-1.d0+1.d-13)then
        print *,' #(E01)# Variable szval given to sz out of range'
        stop
      end if
      if(idim.lt.3)then
        print *,' #(E02)# Incorrect idim or n given to sz'
        stop
      end if
c
c* initialization
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
      iupspn=n/2+mod(n,2)+int(szval+0.001d0)
      icnt=0
      ja=0
      jb=0
      ibpatn=0
c
c* main loop
      do 10 i=1,2**n
        isz=0
        do 20 j=0,n-1
 20     isz=isz+mod(i/2**j,2)
        if(isz.ne.iupspn)go to 10
        icnt=icnt+1
        if(icnt.gt.idim)then
          print *,' #(E02)# Incorrect idim or n given to sz'
          stop
        end if
        ia=iand(i,irght)
        ib=iand(i,ilft)/ihfbit
        if(ib.eq.ibpatn)then
          ja=ja+1
        else
          ibpatn=ib
          ja=1
          jb=icnt-1
        end if
        list1(icnt)=i
        list2(1,ia)=ja
        list2(2,ib)=jb
 10   continue
      if(icnt.eq.idim)return
      print *,' #(E02)# Incorrect idim or n given to sz'
      stop
      end
c
c************* data check of pairs of sites ************
c
      subroutine datack(ipair,ibond,n)
      dimension ipair(ibond*2)
c
      do 10 k=1,ibond
         isite1=ipair(k*2-1)
         isite2=ipair(k*2  )
         if(isite1.le.0.or.isite2.le.0.or.
     &      isite1.gt.n.or.isite2.gt.n)then
            print *,' #(E03)# Incorrect data in ipair'
            print *,'         Location :  ',k*2-1,k*2
            stop
         end if
 10    continue
       end
c
c********* eigenvalues by the bisection method **********
c
      subroutine bisec(alpha,beta,ndim,E,ne,eps)
c
c    alpha  @ diagonal element
c    beta   @ subdiagonal element
c    ndim   @ matrix dimension
c    E      # eigenvalues
c    ne     @ number of eigenvalues to calculate
c    eps    @ limit of error

      implicit real*8 (a-h,o-z)
      dimension alpha(ndim),beta(ndim),E(ne),b2(2000)
c
      if(ndim.gt.2000)then
          print *,' #(E04)# ndim given to bisec exceeds 2000'
          stop
      end if
      if(ne.gt.ndim.or.ne.le.0)then
          print *,' #(E05)# ne given to bisec out of range'
          stop
      end if
c
c*** initial bound
      range=abs(alpha(1))+abs(beta(1))
      do 10 k=2,ndim-1
 10   range=max(range,abs(beta(k-1))+abs(alpha(k))+abs(beta(k)))
      range=max(range,abs(beta(ndim-1))+abs(alpha(ndim)))
      range=-range
c
      b2(1)=0.d0
      do 20 i=2,ndim
 20   b2(i)=beta(i-1)**2
c
      epsabs=abs(range)*eps
      do 30 i=1,ne
 30   E(i)=-range
      b=range
c
c*** bisection method
      do 100 k=1,ne
        a=E(k)
        do 110 j=1,100
          c=(a+b)/2.d0
          if(abs(a-b).lt.epsabs)goto 100
          numneg=0
          g=1.d0
          ipass=0
          do 120 i=1,ndim
            if(ipass.eq.0)then
              g=c-alpha(i)-b2(i)/g
              else if(ipass.eq.1)then
                ipass=2
              else
                g=c-alpha(i)
                ipass=0
            end if
c
            if(ipass.eq.0)then
              if(g.le.0.d0)numneg=numneg+1
              if(abs(g).le.abs(b2(i)*epsabs*eps))ipass=1
            end if
 120      continue
          numneg=ndim-numneg
          if(numneg.lt.k)then
            b=c
          else
            a=c
            do 130 i=k,min(numneg,ne)
 130        E(i)=c
          end if
 110    continue
 100   continue
       end
c
c*** eigenvector of a tridiagonal matrix by inverse iteration ***
c                 for the large/medium routines
c
      subroutine vec12(E,ndim,nvec,di,bl,bu,bv,cm,lex)
c
c    E(4)       @  4 lowest eigenvalues
c    ndim       @  matrix dimension
c    nvec       @  number of vectors to calculate
c    di - lex      working areas
c
      implicit real*8 (a-h,o-z)
      dimension E(4)
      dimension di(ndim),bl(ndim),bu(ndim),bv(ndim),cm(ndim),lex(ndim)
      common /vecdat/alpha(1795),beta(1795),v(1795,5)
c
      do 10 k=1,nvec
c
        do 100 j=1,ndim
          di(j)=E(k)-alpha(j)
          bl(j)=-beta(j)
          bu(j)=-beta(j)
 100    continue
c
c*** LU decomposition
        do 110 j=1,ndim-1
         if(abs(di(j)).gt.abs(bl(j)))then
c--- non pivoting
           lex(j)=0
           if(abs(di(j)).lt.1.d-13)di(j)=1.d-13
           cm(j+1)=bl(j)/di(j)
           di(j+1)=di(j+1)-cm(j+1)*bu(j)
           bv(j)=0.d0
         else
c--- pivoting
           lex(j)=1
           cm(j+1)=di(j)/bl(j)
           di(j)=bl(j)
           s=bu(j)
           bu(j)=di(j+1)
           bv(j)=bu(j+1)
           di(j+1)=s-cm(j+1)*bu(j)
           bu(j+1)= -cm(j+1)*bv(j)
         end if
 110    continue
        if(abs(di(ndim)).lt.1.d-13)di(ndim)=1.d-13
c
c--- initial vector
        do 120 j=1,ndim
 120    v(j,k)=1.d0/(float(j)*5.d0)
c
c*** degeneracy check up
        if(k.eq.1)then
           km=k
        else if(abs(E(k)-E(km)).gt.1.d-13)then
           km=k
        else
           do 130 i=km,k-1
             prd=0.d0
             do 140 j=1,ndim
 140         prd=prd+v(j,i)*v(j,k)
             do 150 j=1,ndim
 150         v(j,k)=v(j,k)-prd*v(j,i)
 130       continue
        end if
c
c*** inverse iteration
        do 160 l=1,k-km+3
          if((l.ne.1).or.(k.ne.km))then
c--- forward substitution
            do 170 j=1,ndim-1
              if(lex(j).eq.0)then
                v(j+1,k)=v(j+1,k)-cm(j+1)*v(j,k)
              else
                s=v(j,k)
                v(j,k)=v(j+1,k)
                v(j+1,k)=s-cm(j+1)*v(j,k)
              end if
 170        continue
          end if
c--- backward substitution
          do 180 j=ndim,1,-1
            s=v(j,k)
            if(j.le.ndim-1)s=s-bu(j)*v(j+1,k)
            if(j.le.ndim-2)s=s-bv(j)*v(j+2,k)
            v(j,k)=s/di(j)
 180      continue
c
c*** normalization
          dnorm=0.d0
          do 190 j=1,ndim
 190      dnorm=dnorm+v(j,k)**2
          if(dnorm.gt.1.d-13)dnorm=1./sqrt(dnorm)
          do 200 j=1,ndim
 200      v(j,k)=v(j,k)*dnorm
 160    continue
c
 10   continue
      end
c
c************* xx correlation function **************
c
      subroutine xcorr(n,idim,npair,nbond,x,sxx,list1,list2)
c
c    n           @ lattice size
c    idim        @ matrix dimension
c    npair       @ pair of sites (k,l) <Sx(k)Sx(l)>
c    nbond       @ number of bonds to be investigated
c    x           @ eigenvetor
c    sxx         # xx correlation function
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8 (a-h,o-z)
      dimension list1(idim),list2(2,0:2**24)
      dimension x(idim)
      dimension npair(nbond*2),sxx(nbond)
c
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
c
      do 10 k=1,nbond
        i1=npair(k*2-1)-1
        i2=npair(k*2  )-1
        if(i1.lt.0.or.i1.ge.n.or.i2.lt.0.or.i2.ge.n.or.i1.eq.i2)then
          print *,' #(W01)# Wrong site number given to xcorr'
          return
        end if
        corr=0.d0
        is=2**i1+2**i2
!$OMP parallel do 
!$OMP& private(j,ibit,iexchg,ia,ib),
!$OMP& shared(idim,list1,is,irght,ilft,ihfbit,x,list2),
!$OMP& reduction(+:corr)
        do 20 j=1,idim
          ibit=iand(list1(j),is)
          if(ibit.eq.0.or.ibit.eq.is)goto 20
          iexchg=ieor(list1(j),is)
          ia=iand(iexchg,irght)
          ib=iand(iexchg,ilft)/ihfbit
          corr=corr+x(j)*x(list2(1,ia)+list2(2,ib))
 20     continue
!$OMP end parallel do 
        sxx(k)=corr/4.d0
 10   continue
      return
      end
c
c************* zz correlation function **************
c
      subroutine zcorr(n,idim,npair,nbond,x,szz,list1)
c
c    n           @ lattice size
c    idim        @ matrix dimension
c    npair       @ pair of sites (k,l) <Sz(k)Sz(l)>
c    nbond       @ number of bonds to be investigated
c    x           @ eigenvetor
c    szz         # zz correlation function
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8 (a-h,o-z)
      dimension list1(idim)
      dimension x(idim)
      dimension npair(nbond*2),szz(nbond)
c
      do 10 k=1,nbond
        i1=npair(k*2-1)-1
        i2=npair(k*2  )-1
        if(i1.lt.0.or.i1.ge.n.or.i2.lt.0.or.i2.ge.n.or.i1.eq.i2)then
          print *,' #(W02)# Wrong site number given to zcorr'
          return
        end if
        corr=0.d0
        is=2**i1+2**i2
!$OMP parallel do private(j,ibit,factor),shared(idim,list1,is,x),
!$OMP& reduction(+:corr) 
        do 20 j=1,idim
          ibit=iand(list1(j),is)
          if(ibit.eq.0.or.ibit.eq.is)then
            factor=1.d0
          else
            factor=-1.d0
          end if
          corr=corr+factor*x(j)**2
 20     continue
!$OMP end parallel do
        szz(k)=corr/4.d0
 10   continue
      return
      end
c
c******* Orthogonalization of the eigenvectors ************
c
      subroutine orthg(idim,ideclr,ev,norm,idgn,numvec)
c
c   idim    @  matrix dimension
c   ideclr  @  declared array size in the main program
c   ev      @# vectors to be orthogonalized / orthogonalized vectors
c   norm(j) #  norm of the j-th vector returned
c   idgn    #  degree of degenearcy
c   numvec  @  number of vectors to be checked
c
      implicit real*8(a-h,o-z)
      dimension ev(ideclr,numvec),norm(numvec)
      if(numvec.le.1)then
         print *,' #(W03)# Number of vectors is less than 2 in orthg'
         return
      end if
      do 10 i=1,numvec
        dnorm=0.0d0
!$OMP parallel do private(j),shared(idim,ev,i),reduction(+:dnorm)
        do 20 j=1,idim
 20     dnorm=dnorm+ev(j,i)**2
!$OMP end parallel do
        if(dnorm.lt.1.0d-20)then
           print *,' #(W04)# Null vector given to orthg. Location is',i
           return
        end if
        dnorm=1.0d0/sqrt(dnorm)
!$OMP parallel do private(j),shared(idim,ev,dnorm)
        do 25 j=1,idim
 25     ev(j,i)=ev(j,i)*dnorm
!$OMP end parallel do
 10   continue
      idgn=numvec
      norm(1)=1
c*** orthogonalization
      do 30 i=2,numvec
       norm(i)=1
       do 40 j=1,i-1
         prjct=0.0d0
!$OMP parallel do private(l),shared(idim,ev,j),reduction(+:prjct)
         do 50 l=1,idim
 50      prjct=prjct+ev(l,i)*ev(l,j)
!$OMP end parallel do
!$OMP parallel do private(l),shared(idim,ev,j,prjct)
         do 60 l=1,idim
 60      ev(l,i)=ev(l,i)-prjct*ev(l,j)
!$OMP end parallel do
 40    continue
       vnorm=0.0d0
!$OMP parallel do private(l),shared(idim,ev),reduction(+:vnorm)
       do 70 l=1,idim
 70    vnorm=vnorm+ev(l,i)**2
!$OMP end parallel do
       if(vnorm.gt.1.0d-15)then
         vnorm=1.0d0/sqrt(vnorm)
!$OMP parallel do private(l),shared(idim,ev,vnorm)
         do 80 l=1,idim
 80      ev(l,i)=ev(l,i)*vnorm
!$OMP end parallel do
        else
!$OMP parallel do private(l),shared(idim,ev)
         do 90 l=1,idim
 90      ev(l,i)=0.0d0
!$OMP end parallel do
         idgn=idgn-1
         norm(i)=0
       end if
 30   continue
c*** check orthogonality
      do 100 i=2,numvec
       do 100 j=1,i-1
       prd=0.0d0
!$OMP parallel do private(l),shared(idim,ev),reduction(+:prd)
       do 110 l=1,idim
 110   prd=prd+ev(l,i)*ev(l,j)
!$OMP end parallel do
       if(abs(prd).lt.1.0d-10)go to 100
         print 200,i,j
 200     format(' #(W05)# Non-orthogonal vectors at',2i4)
         print 210,prd
 210     format('         Overlap : ',d14.7)
         print *,'         Unsuccessful orthogonalization'
         return
 100  continue
      return
      end
c
c***********************************************************
c*                                                         *
c*                   TITPACK Ver. 2                        *
c*                                                         *
c*                   February, 1991                        *
c*                                                         *
c*         Copyright (C) Hidetoshi Nishimori               *
c*                                                         *
c***********************************************************
c
c============  SUBROUTINES / LARGE MATRICES =================
c
c*** variables marked @ should be given from the main program
c*** variables marked # are evaluated and returned
c*** The following variables are common to all routines
c    n         @ lattice size
c    idim      @ matrix dimension
c    ipair     @ pairs of sites connected by bonds
c    bondwt    @ exchange interaction of each bond Jxy
c    zrtio     @ ratio of Jz to Jxy
c    ibond     @ number of bonds
c
c************ eigenvalues by the Lanczos method **************
c         --- dummy routine for simple working area
c
      subroutine lnc1(n,idim,ipair,bondwt,zrtio,ibond,
     &                nvec,iv,E,itr,wk,ideclr,list1,list2)
c
c    nvec      @ number of eigenvectors to calculate in lncvec
c    iv        @ location of the nonzero element of the initial vector
c    E         # eigenvalues
c    itr     # number of iterations required for convergence
c    wk       working areas
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension E(4)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**24)
      dimension wk(ideclr,2)
c
      if(iv.le.0.or.iv.gt.idim)then
          print *,' #(E06)# Incorrect iv given to lnc1'
          return
      end if
      if(nvec.lt.0.or.nvec.gt.4)then
          print *,' #(W06)# Wrong value given to nvec in lnc1'
          print *,'         Only the eigenvalues are calculated'
          nvec=0
      end if
c
      call lnc1z(n,idim,ipair,bondwt,zrtio,ibond,
     &     nvec,iv,E,itr,wk(1,1),wk(1,2),list1,list2)
      end
c
c************ eigenvalues by the Lanczos method
c
      subroutine lnc1z(n,idim,ipair,bondwt,zrtio,ibond,
     &                nvec,iv,E,itr,v1,v0,list1,list2)
c
      implicit real*8(a-h,o-z)
      dimension E(4)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**24)
      dimension v0(idim),v1(idim)
      dimension wk(1795,5),iwk(1795)
      common /vecdat/alpha(1795),beta(1795),coef(1795,5)
c
c*** initialization
!$OMP parallel do private(i),shared(idim,v0,v1)
      do 10 i=1,idim
          v0(i)=0.0d0
          v1(i)=0.0d0
 10   continue
!$OMP end parallel do
c
      v1(iv)=1.0d0
c
      call datack(ipair,ibond,n)
c
c*** alpha(1) and beta(1)
      call mltply(n,idim,ipair,bondwt,zrtio,ibond,
     &            v1,v0,prdct,list1,list2)
      alpha1=prdct
      alpha(1)=alpha1
      beta1=0.d0
!           write(*,*) 1,"-st step: ",alpha1,beta1
!$OMP parallel do private(i),shared(idim,v0,alpha1,v1),
!$OMP& reduction(+:beta1)
      do 50 i=1,idim
 50   beta1=beta1+(v0(i)-alpha1*v1(i))**2
!$OMP end parallel do
      beta1=sqrt(beta1)
      beta(1)=beta1
c
C*** iteration
      do 100 i=2,1795
!           write(*,*) i,"-th step: ",alpha1,beta1
!$OMP parallel do private(j,temp1,temp2),
!$OMP& shared(idim,v1,v0,alpha1,beta1)
          do 110 j=1,idim
            temp1=v1(j)
            temp2=(v0(j)-alpha1*v1(j))/beta1
            v0(j)=-beta1*temp1
            v1(j)=temp2
 110      continue
!$OMP end parallel do
c
          call mltply(n,idim,ipair,bondwt,zrtio,ibond,
     &                v1,v0,prdct,list1,list2)
          alpha1=prdct
          alpha(i)=alpha1
          beta1=0.d0
!$OMP parallel do private(j),shared(idim,v0,alpha1,v1),
!$OMP& reduction(+:beta1)
          do 120 j=1,idim
 120      beta1=beta1+(v0(j)-alpha1*v1(j))**2
!$OMP end parallel do
          beta1=sqrt(beta1)
          beta(i)=beta1
          if(beta(i).lt.0.5d-30)then
            print *,' #(E07)# Tridiagonalization unsuccessful in lnc1'
            print *,'         Beta(i) is too small at i=  ',i
            stop
          end if
c
c*** convergence check
          if(i.gt.20.and.mod(i,5).eq.0)then
             call bisec(alpha,beta,i,E,4,eps)
           write(*,*) i,"-th energies: ",e(1),e(2)
             if(abs((ebefor-E(2))/E(2)).lt.1.0d-13)then
                if(nvec.gt.0)call vec12(E,i,nvec,wk(1,1),wk(1,2),
     &                            wk(1,3),wk(1,4),wk(1,5),iwk)
                itr=i
                return
             end if
             ebefor=E(2)
          end if
          if(i.eq.20)then
             eps=1.d-10
             call bisec(alpha,beta,20,E,4,eps)
             ebefor=E(2)
          end if
 100  continue
c
      print *,' #(W07)# lnc1 did not converge within 1795 steps'
      itr=1795
      return
      end
c
c************ eigenvector by the Lanczos method *************
c
      subroutine lncv1(n,idim,ipair,bondwt,zrtio,ibond,
     &                 nvec,iv,x,itr,wk,ideclr,list1,list2)
c
c    nvec      @ number of eigenvectors to be calculated
c    iv        @ location of the nonzero element of the initial vector
c    x         # eigenvector
c    itr     @ number of interations for convergence
c    wk         working area
c    ideclr    @ declared array size in the main program
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension x(ideclr,nvec)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**24)
      dimension wk(ideclr,2)
c
      if(nvec.le.0.or.nvec.gt.4)then
          print *,'#(W08)# nvec given to lncv1 out of range'
          return
      end if
      call lncv1z(n,idim,ipair,bondwt,zrtio,ibond,
     &       nvec,iv,x,ideclr,itr,wk(1,1),wk(1,2),list1,list2)
      end
c
c************ eigenvector by the Lanczos method
c
      subroutine lncv1z(n,idim,ipair,bondwt,zrtio,ibond,
     &                 nvec,iv,x,ideclr,itr,v1,v0,list1,list2)
c
      implicit real*8(a-h,o-z)
      dimension x(ideclr,nvec)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**24)
      dimension v0(idim),v1(idim)
      common /vecdat/alpha(1795),beta(1795),coef(1795,5)
c
c*** initialization
!$OMP parallel do private(i),shared(idim,v0,v1)
      do 10 i=1,idim
          v0(i)=0.0d0
          v1(i)=0.0d0
 10   continue
!$OMP end parallel do
      do 15 k=1,nvec
!$OMP parallel do private(i),shared(idim,x,k)
      do 14 i=1,idim
 14   x(i,k)=0.d0
!$OMP end parallel do
 15   continue
c
      v1(iv)=1.0d0
      do 20 k=1,nvec
 20   x(iv,k)=coef(1,k)
c
c*** alpha(1) and beta(1)
      call mltply(n,idim,ipair,bondwt,zrtio,ibond,
     &            v1,v0,prdct,list1,list2)

      alpha1=alpha(1)
      beta1=beta(1)
!           write(*,*) 1,"-st step: ",alpha1,beta1
      do 40 k=1,nvec
!$OMP parallel do private(j),shared(idim,x,coef,v0,alpha1,v1,beta1)
      do 39 j=1,idim
 39   x(j,k)=x(j,k)+coef(2,k)*(v0(j)-alpha1*v1(j))/beta1
!$OMP end parallel do
 40   continue
c
c*** iteration
      do 100 i=2,itr-1
!           write(*,*) i,"-th step: ",alpha1,beta1
!$OMP parallel do private(j,temp1,temp2),shared(idim,v1,v0,alpha1,beta1)
          do 110 j=1,idim
            temp1=v1(j)
            temp2=(v0(j)-alpha1*v1(j))/beta1
            v0(j)=-beta1*temp1
            v1(j)=temp2
 110      continue
!$OMP end parallel do
c
          call mltply(n,idim,ipair,bondwt,zrtio,ibond,
     &                v1,v0,prdct,list1,list2)
          alpha1=alpha(i)
          beta1=beta(i)
          do 130 k=1,nvec
!$OMP parallel do private(j),shared(idim,x,k,coef,v0,alpha1,v1,beta1)
          do 129 j=1,idim
 129      x(j,k)=x(j,k)+coef(i+1,k)*(v0(j)-alpha1*v1(j))/beta1
!$OMP end parallel do
 130      continue
 100  continue
c
c*** normalization
      do 200 k=1,nvec
          dnorm=0.d0
!$OMP parallel do private(j),shared(idim,x),reduction(+:dnorm)
          do 210 j=1,idim
 210      dnorm=dnorm+x(j,k)**2
!$OMP end parallel do
          dnorm=sqrt(dnorm)
!$OMP parallel do private(j),shared(idim,x,dnorm)
          do 220 j=1,idim
 220      x(j,k)=x(j,k)/dnorm
!$OMP end parallel do
 200  continue
      return
      end
c
c*************** matrix multiplication *****************
c
      subroutine mltply(n,idim,ipair,bondwt,zrtio,ibond,
     &                  v1,v0,prdct,list1,list2)
c
c    v1          @  input vector
c    v0          @# output vector : H*v1+v0(input)
c    prdct       #  <v1*H*v1>
c    list1,list2 @  spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**24)
      dimension v0(idim),v1(idim)
c
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
c
      prdct=0.d0
      do 10 k=1,ibond
          isite1=ipair(k*2-1)-1
          isite2=ipair(k*2  )-1
          is1=2**isite1
          is2=2**isite2
          is=is1+is2
          wght=bondwt(k)*zrtio(k)*0.5d0
!$OMP parallel do private(j,ibit,factor,offdg,iexchg,ia,ib,temp),
!$OMP& shared(idim,list1,is,v0,wght,v1,irght,ilft,ihfbit,list2,bondwt,k),
!$OMP& reduction(+:prdct)
          do 20 j=1,idim
            ibit=iand(list1(j),is)
            if(ibit.eq.0.or.ibit.eq.is)then
               v0(j)=v0(j)-wght*v1(j)
               factor=1.d0
               offdg=0.d0
            else
               v0(j)=v0(j)+wght*v1(j)
               iexchg=ieor(list1(j),is)
               ia=iand(iexchg,irght)
               ib=iand(iexchg,ilft)/ihfbit
               temp=v1(list2(1,ia)+list2(2,ib))*bondwt(k)
               v0(j)=v0(j)-temp
               factor=-1.d0
               offdg=-temp*v1(j)
            end if
            prdct=prdct-factor*wght*v1(j)**2+offdg
 20   continue
!$OMP end parallel do
 10   continue
      end
c
c*************** check of the eigenvector and eigenvalue ************
c
      subroutine check1(n,idim,ipair,bondwt,zrtio,ibond,
     &                  x,v,Hexpec,list1,list2)
c
c    x           @ eigenvector to be checked
c    v           # H*x
c    Hexpec      # <x*H*x>
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension x(idim),v(idim)
      dimension list1(idim),list2(2,0:2**24)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
c
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
c
      dnorm=0.d0
!$OMP parallel do private(i),shared(idim,x),reduction(+:dnorm)
      do 5 i=1,idim
 5    dnorm=dnorm+x(i)**2
!$OMP end parallel do
      if(dnorm.lt.1.d-30)then
         print *,' #(W09)# Null vector given to check1'
         return
       end if
!$OMP parallel do private(i),shared(idim,v)
      do 10 i=1,idim
 10   v(i)=0.0d0
!$OMP end parallel do
      do 20 k=1,ibond
          isite1=ipair(k*2-1)-1
          isite2=ipair(k*2  )-1
          is1=2**isite1
          is2=2**isite2
          is=is1+is2
          wght=bondwt(k)
          diag=wght*0.5d0*zrtio(k)
!$OMP parallel do private(i,ibit,iexchg,ia,ib),
!$OMP& shared(idim,list1,is,v,diag,x,irght,ilft,ihfbit,list2,wght)
          do 2000 i=1,idim
            ibit=iand(list1(i),is)
             if(ibit.eq.0.or.ibit.eq.is)then
                v(i)=v(i)-diag*x(i)
               else
                v(i)=v(i)+diag*x(i)
                iexchg=ieor(list1(i),is)
                ia=iand(iexchg,irght)
                ib=iand(iexchg,ilft)/ihfbit
                v(i)=v(i)-
     &               x(list2(1,ia)+list2(2,ib))*wght
             end if
 2000     continue
!$OMP end parallel do
 20   continue
c
      prd=0.0d0
!$OMP parallel do private(i),shared(idim,v,x),reduction(+:prd)
      do 30 i=1,idim
 30   prd=prd+v(i)*x(i)
!$OMP end parallel do
      Hexpec=prd
      print *
      print 200
 200  format(' ---------------------------- Information from check1')
      print 210,prd
 210  format(' <x*H*x> =',1pd16.8)
      print 220
 220  format(' H*x(j)/x(j) (j=min(idim/3,13),idim,max(1,idim/20))')
      print 230,(v(i)/x(i),i=min(idim/3,13),idim,max(1,idim/20))
 230  format(4d18.9)
      print 240
 240  format(' ---------------------------------------------------')
      return
      end
c
c****************** inverse iteration ************************
c      --- dummy routine for simple working area
c
      subroutine inv1(n,idim,ipair,bondwt,zrtio,ibond,
     &                  Eig,iv,x,wk,ideclr,list1,list2)
c
c    Eig         @ eigenvalue
c    iv          @ location of the nonzero element of the initial vector
c    x           # eigenvector
c    wk            working area
c    ideclr      @ declared array size in the main program
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**24)
      dimension x(idim)
      dimension wk(ideclr,4)
c
      call inv1z(n,idim,ipair,bondwt,zrtio,ibond,
     & Eig,iv,x,wk(1,1),wk(1,2),wk(1,3),wk(1,4),list1,list2)
      end
c
c****************** inverse iteration
c
      subroutine inv1z(n,idim,ipair,bondwt,zrtio,ibond,
     &                  Eig,iv,x,b,p,r,y,list1,list2)
c
c    b            working area for the rhs of (H-E(approx))*x=b
c    p,r,y        working area used in the routine cg1
c


      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**24)
      dimension b(idim),x(idim),r(idim),y(idim),p(idim)
c
      do 10 i=1,idim
 10   b(i)=0.0d0
      b(iv)=1.0d0
      do 20 itr=1,20
        call cg1(n,idim,ipair,bondwt,zrtio,ibond,
     &          Eig,x,b,p,r,y,iterat,list1,list2)
        if(iterat.gt.idim)then
          xnorm=0.0d0
          do 22 i=1,idim
 22       xnorm=xnorm+x(i)**2
          xnorm=sqrt(xnorm)
          do 24 i=1,idim
 24       x(i)=x(i)/xnorm
          print *,' #(W10)# Iterat in cg1 exceeds idim or 500'
          print *,'         Approximate eigenvector returned'
          print *,'         Itration number in inv1 is',itr
          return
        end if
        xnorm=0.0d0
        do 30 i=1,idim
 30     xnorm=xnorm+x(i)**2
        xnorm=sqrt(xnorm)
        do 40 i=1,idim
 40     x(i)=x(i)/xnorm
        xb=0.0d0
        do 50 i=1,idim
 50     xb=xb+x(i)*b(i)
        if(abs(abs(xb)-1.0d0).lt.1.0d-12)then
c          print 100,itr
 100      format('       number of iterations in inv1 :',i5)
          return
        end if
        do 60 i=1,idim
 60     b(i)=x(i)
 20   continue
      print *,' #(W11)# inv1 did not converge'
      return
      end
c
c************** solution of linear equations -- cg method ************
c
      subroutine cg1(n,idim,ipair,bondwt,zrtio,ibond,
     &              Eig,x,b,p,r,y,itr,list1,list2)
c
c    Eig         @ eigenvalue
c    x           # eigenvector
c    b             working area for the rhs of (H-E(approx))*x=b
c    p,r,y         working area used in the routine cg
c    itr         # number of iterations required for convergence
c    list1,list2 @ spin configurations generated in 'sz'
c
      implicit real*8(a-h,o-z)
      dimension ipair(ibond*2),bondwt(ibond),zrtio(ibond)
      dimension list1(idim),list2(2,0:2**24)
      dimension b(idim),x(idim),r(idim),y(idim),p(idim)
c
      ihfbit=2**((n+1)/2)
      irght=2**((n+1)/2)-1
      ilft=ieor(2**n-1,irght)
c
c*** initialization
      bnorm=0.0d0
      do 10 i=1,idim
          bnorm=bnorm+b(i)**2
          r(i)=b(i)
          p(i)=b(i)
          x(i)=0.0d0
 10   continue
c
c*** iteration
      do 20 itr=1,min(500,idim)
        do 30 i=1,idim
 30     y(i)=0.0d0
        do 40 k=1,ibond
          isite1=ipair(k*2-1)-1
          isite2=ipair(k*2  )-1
          is1=2**isite1
          is2=2**isite2
          is=is1+is2
          eperbd=Eig/float(ibond)
          wght=bondwt(k)
          diag1= wght*0.5d0*zrtio(k)+eperbd
          diag2=-wght*0.5d0*zrtio(k)+eperbd
          do 40 i=1,idim
            ibit=iand(list1(i),is)
            if(ibit.eq.0.or.ibit.eq.is)then
                y(i)=y(i)-diag1*p(i)
            else
                iexchg=ieor(list1(i),is)
                ia=iand(iexchg,irght)
                ib=iand(iexchg,ilft)/ihfbit
                y(i)=y(i)-diag2*p(i)-p(list2(1,ia)+list2(2,ib))*wght
            end if
 40     continue
        rp=0.0d0
        yp=0.0d0
        do 50 i=1,idim
          rp=rp+r(i)*p(i)
          yp=yp+y(i)*p(i)
 50     continue
        alpha=rp/yp
        rnorm=0.0d0
        do 60 i=1,idim
          x(i)=x(i)+alpha*p(i)
          rnorm=rnorm+r(i)**2
 60     continue
        rnorm2=0.0d0
        do 70 i=1,idim
          r(i)=r(i)-alpha*y(i)
          rnorm2=rnorm2+r(i)**2
 70     continue
        beta=rnorm2/rnorm
        do 90 i=1,idim
 90     p(i)=r(i)+beta*p(i)
        if(mod(itr,5).ne.0)go to 20
        if(sqrt(rnorm2).lt.1.0d-9*sqrt(bnorm))then
c          print 150,itr
 150      format('       number of iterations in cg1     :',i5)
          return
        end if
 20   continue
c     print *,' #(Wxx)# cg1 did not converge'
      return
      end
c************* z correlation function **************
c
      subroutine zexpect(n,idim,isite,x,szexp,list1)
c
c    n           @ lattice size
c    idim        @ matrix dimension
c    isite       @ site k for <Sz(k)>
c    x           @ eigenvetor
c    szexp       # expectation value of Sz(k)
c    list1,list2 @ spin configurations generated in 'sz'
c
c      implicit real*8 (a-h,o-z)
      implicit none
c input
      integer n
      integer idim
      integer list1
      dimension list1(idim)
      real*8 x
      dimension x(idim)
      integer isite
c output
      real*8 szexp
c local
      integer j
      real*8 corr
      integer i1,is,ibit
      real*8 factor
c
        i1=isite-1
        if((i1.lt.0).or.(i1.ge.n))then
          print *,' #(W02)# Wrong site number given to zexp'
          return
        end if
        corr=0.d0
        is=2**i1
!$OMP parallel do
!$OMP& private(j,ibit,factor),
!$OMP& shared(idim,list1,is,x),
!$OMP& reduction(+:corr)
        do 20 j=1,idim
          ibit=iand(list1(j),is)
          if(ibit.eq.is)then
            factor=1.d0
          else
            factor=-1.d0
          end if
          corr=corr+factor*x(j)**2
 20                   continue
!$OMP end parallel do
c
        szexp=corr/2.d0
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine name_sub
c      implicit none
cc
cc  input parameters
cc output parameters
cc local variables
cc begin
c      return
c      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc