!*********************************
!* Sample program -- S=1/2 chain *
!* Date:2007.5.17(Thu)           *
!*********************************

!**********
 module com
!**********
 integer,parameter::NSITE=8,NPAIR=8
 real*8,parameter::CEPS=1.d-13
 end module com

!*************
 program main
 use com
 implicit none
!*************
 integer::dim,ne,nvec,ns,np,i,pair(0:NPAIR-1,2)
 real*8::dummy(1),eps
 real*8,dimension(:,:),allocatable::element,wk
 real*8,dimension(:),allocatable::ev
 integer,dimension(:),allocatable::iwk

   ns=NSITE
   dim=2**ns
   np=NPAIR
   call set_pair(ns,np,pair)

   allocate(element(dim,dim))
   allocate(wk(dim,8))
   allocate(ev(dim))
   allocate(iwk(dim))
   
   call set_element(dim,np,pair,element)
   ne=dim
   nvec=0
   eps=CEPS
   call diag(element,dim,dim,ev,dummy,ne,nvec,eps,wk,iwk)

   write(*,'("ground state energy:",E23.15)') ev(1)*2.D0

   open(1,FILE="sample.dat")
   write(1,'("ns=",I3)') ns
   write(1,'(I8)') dim
   do i=1,dim
      write(1,'(E23.15)') ev(i)
   enddo
   close(1)

   open(3,FILE="eightsite.dat")
   write(3,'("ground state energy:",E23.15)') ev(1)*2.D0
   close(3)

   deallocate(element)
   deallocate(wk)
   deallocate(ev)
   deallocate(iwk)

 end program main

!****************************************
 subroutine set_pair(ns,np,pair) !07.5.17
 implicit none
 integer,intent(in)::ns,np
 integer,intent(out)::pair(0:np-1,2)
!****************************************
 integer::i

   do i=0,ns-2
      pair(i,:)=[i,i+1]
   enddo
   pair(ns-1,:)=[ns-1,0]

 end subroutine set_pair

!****************************************************
 subroutine set_element(dim,np,pair,element) !07.5.17
 implicit none
 integer,intent(in)::dim,np,pair(0:np-1,2)
 real*8,intent(out)::element(dim,dim)
 integer disp,indisp
!****************************************************
 integer::k1,k2,l,i_a,i_b,s_a,s_b,mask

   element=0.D0
   do k1=0,dim-1
      do l=0,np-1
         i_a=pair(l,1)
         s_a=ibits(k1,i_a,1)
         i_b=pair(l,2)
         s_b=ibits(k1,i_b,1)
         element(k1+1,k1+1)=element(k1+1,k1+1)+dble((2*s_a-1)*(2*s_b-1))*0.25d0
         if(s_a/=s_b) then
            mask=ibset(ibset(0,i_a),i_b)
            k2=ieor(k1,mask)
            !排他的論理和
            element(k2+1,k1+1)=element(k2+1,k1+1)+0.5d0
         endif
      enddo
   enddo
   open(2,file = "hamirutonian.txt")
   do disp = 1,dim
      do indisp = 1,dim
            write(2,'(F10.5)',advance="no")element(disp,indisp)
            write(2, '(A)', advance="no") " "   ! スペースを追加
      end do 
      write(2, '(A)')""
   end do
   close(2)

 end subroutine set_element

!************ eigenvalues of a small matrix *************                       
!                                                                               
subroutine diag(elemnt,ideclr,idim,E,v,ne,nvec,eps,wk,iwk)                
!                                                                               
!    elemnt      @ matrix elements                                              
!    ideclr      @ declared array size in the main program                      
!    idim        @ matrix dimension                                             
!    E           # eigenvalues                                                  
!    v           # eigenvector                                                  
!    ne          @ number of eigenvalues to calculate                           
!    nvec        @ number of eigenvectors to calculate                          
!    eps         @ limit of error                                               
!    wk,iwk        working areas                                                
!                                                                               
      implicit real*8(a-h,o-z)                                                  
      dimension E(ne),v(ideclr,nvec),elemnt(ideclr,idim)                        
      dimension wk(ideclr,8),iwk(ideclr)                                        
!                                                                               
      if(nvec.lt.0.or.nvec.gt.ne)then                                           
           print *,' #(E10)# nvec given to diag out of range'                   
           stop                                                                 
      end if                                                                    
!                                                                               
      call hshldr(elemnt,ideclr,idim,wk(1,1),wk(1,2),wk(1,3),wk(1,4),&
     &            wk(1,5),wk(1,6))                                              
      call bisec(wk(1,1),wk(1,2),idim,E,ne,eps)                                 
      if(nvec.eq.0)return                                                       
      call vec3(E,elemnt,ideclr,idim,ne,nvec,wk(1,4),wk(1,5),wk(1,6),&
     &          wk(1,7),wk(1,8),iwk,wk(1,1),wk(1,2),wk(1,3),v)                  
      end                                                                       
!                                                                               
!************ Householder tridiagonalization *************                      
!                                                                               
      subroutine hshldr(elemnt,ideclr,idim,alpha,beta,c,w,p,q)                  
!                                                                               
!    elemnt     @ matrix                                                        
!    ideclr     @ declared array size in the main program                       
!    idim       @ matrix dimension                                              
!    alpha      # diagonal element                                              
!    beta       # subdiagonal element                                           
!    c          # normalization factor for eigenvector calculation              
!    w,p,q        working areas                                                 
!                                                                               
      implicit real*8(a-h,o-z)                                                  
      dimension elemnt(ideclr,idim),alpha(idim),beta(idim),c(idim)              
      dimension w(idim),p(idim),q(idim)                                         
!                                                                               
      do 10 k=1,idim-2                                                          
        s=0.d0                                                                  
        do 20 i=k+1,idim                                                        
 20     s=s+elemnt(i,k)**2                                                      
        s=sqrt(s)                                                               
        if(elemnt(k+1,k).lt.0.d0)s=-s                                           
!                                                                               
        alpha(k)=elemnt(k,k)                                                    
        beta(k)=-s                                                              
        c(k)=0.0d0                                                              
        if(s**2.lt.1.d-26)goto 10                                               
!                                                                               
        c(k)=1.d0/(s**2+elemnt(k+1,k)*s)                                        
        w(k+1)=elemnt(k+1,k)+s                                                  
        do 30 i=k+2,idim                                                        
 30     w(i)=elemnt(i,k)                                                        
        elemnt(k+1,k)=w(k+1)                                                    
!                                                                               
        do 40 i=k+1,idim                                                        
          t=0.d0                                                                
          do 50 j=k+1,i                                                         
 50       t=t+elemnt(i,j)*w(j)                                                  
          do 55 j=i+1,idim                                                      
 55       t=t+elemnt(j,i)*w(j)                                                  
          p(i)=c(k)*t                                                           
 40     continue                                                                
!                                                                               
        t=0.d0                                                                  
        do 60 i=k+1,idim                                                        
 60     t=t+p(i)*w(i)                                                           
        s=0.5d0*c(k)*t                                                          
        do 70 i=k+1,idim                                                        
 70     q(i)=p(i)-s*w(i)                                                        
!                                                                               
        do 80 j=k+1,idim                                                        
        do 80 i=j,idim                                                          
 80     elemnt(i,j)=elemnt(i,j)-w(i)*q(j)-q(i)*w(j)                             
!                                                                               
 10   continue                                                                  
!                                                                               
      alpha(idim-1)=elemnt(idim-1,idim-1)                                       
      alpha(idim)=elemnt(idim,idim)                                             
      beta(idim-1)=elemnt(idim,idim-1)                                          
      return                                                                    
      end                                                                       
!                                                                               
!**** eigenvector of a tridiagonal matrix by inverse iteration ****             
!                                                                               
      subroutine vec3(E,elemnt,ideclr,idim,ne,nvec,di,bl,bu,bv,cm,lex,&
     &                alpha,beta,c,v)                                           
!                                                                               
!    E        @ eigenvalues                                                     
!    elemnt   @ matrix elements                                                 
!    ideclr   @ declared array size in the main program                         
!    idim     @ matrix dimension                                                
!    ne       @ number of eigenvalues                                           
!    nvec     @ number of eigenvectors                                          
!    di - lex   working areas                                                   
!    alpha    @ diagonal element in the tridiagonal form                        
!    beta     @ subdiagonal element in the tridiagonal form                     
!    c        @ normalization factor given from hshldr                          
!    v        # eigenvectors                                                    
!                                                                               
      implicit real*8 (a-h,o-z)                                                 
      dimension E(ne),elemnt(ideclr,idim)                                       
      dimension di(idim),bl(idim),bu(idim),bv(idim),cm(idim),lex(idim)          
      dimension alpha(idim),beta(idim),c(idim),v(ideclr,nvec)                   
!                                                                               
      if(nvec.gt.ne)then                                                        
          print *,' #(E11)# nvec given to vec3 out of range'                    
          stop                                                                  
      end if                                                                    
      do 10 k=1,nvec                                                            
!                                                                               
        do 100 j=1,idim                                                         
          di(j)=E(k)-alpha(j)                                                   
          bl(j)=-beta(j)                                                        
          bu(j)=-beta(j)                                                        
 100    continue                                                                
!                                                                               
!*** LU decomposition                                                           
        do 110 j=1,idim-1                                                       
         if(abs(di(j)).gt.abs(bl(j)))then                                       
!--- non pivoting                                                               
           lex(j)=0                                                             
           if(abs(di(j)).lt.1.d-13)di(j)=1.d-13                                 
           cm(j+1)=bl(j)/di(j)                                                  
           di(j+1)=di(j+1)-cm(j+1)*bu(j)                                        
           bv(j)=0.d0                                                           
         else                                                                   
!--- pivoting                                                                   
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
        if(abs(di(idim)).lt.1.d-13)di(idim)=1.d-13                            
!                                                                               
!--- initial vector                                                             
        do 120 j=1,idim                                                         
 120    v(j,k)=1.d0/dble(j*5.d0)                                               
!                                                                               
!*** degeneracy check up                                                        
        if(k.eq.1)then                                                          
           km=k                                                                 
        else if(abs(E(k)-E(km)).gt.1.d-13)then                                  
           km=k                                                                 
        else                                                                    
           do 130 i=km,k-1                                                      
             prd=0.d0                                                           
             do 140 j=1,idim                                                    
 140         prd=prd+v(j,i)*v(j,k)                                              
             do 150 j=1,idim                                                    
 150         v(j,k)=v(j,k)-prd*v(j,i)                                           
 130       continue                                                             
        end if                                                                  
!                                                                               
!*** inverse iteration                                                          
        do 160 l=1,k-km+3                                                       
          if((l.ne.1).or.(k.ne.km))then                                         
!--- forward substitution                                                       
            do 170 j=1,idim-1                                                   
              if(lex(j).eq.0)then                                               
                v(j+1,k)=v(j+1,k)-cm(j+1)*v(j,k)                                
              else                                                              
                s=v(j,k)                                                        
                v(j,k)=v(j+1,k)                                                 
                v(j+1,k)=s-cm(j+1)*v(j,k)                                       
              end if                                                            
 170        continue                                                            
          end if                                                                
!-- backward substitution                                                      
          do 180 j=idim,1,-1                                                    
            s=v(j,k)                                                            
            if(j.le.idim-1)s=s-bu(j)*v(j+1,k)                                   
            if(j.le.idim-2)s=s-bv(j)*v(j+2,k)                                   
            v(j,k)=s/di(j)                                                      
 180      continue                                                              
!                                                                               
!*** normalization                                                              
          dnorm=0.d0                                                            
          do 190 j=1,idim                                                       
 190      dnorm=dnorm+v(j,k)**2                                                 
          if(dnorm.gt.1.d-13)dnorm=1./sqrt(dnorm)                               
          do 200 j=1,idim                                                       
 200      v(j,k)=v(j,k)*dnorm                                                   
 160    continue                                                                
!                                                                               
 10   continue                                                                  
!                                                                               
!*** back transformation to the original representation                         
      do 210 k=1,nvec                                                           
!                                                                               
        do 220 i=idim-2,1,-1                                                    
          prd=0.d0                                                              
          do 230 j=i+1,idim                                                     
 230      prd=prd+elemnt(j,i)*v(j,k)                                            
          s=prd*c(i)                                                            
          do 240 j=i+1,idim                                                     
 240      v(j,k)=v(j,k)-s*elemnt(j,i)                                           
 220    continue                                                                
 210  continue                                                                  
!                                                                               
!*** orthogonalization for degenerate case                                      
      km=1                                                                      
      do 250 k=2,nvec                                                           
        if(abs(E(k)-E(km)).ge.1.0d-13)then                                      
           km=k                                                                 
        else                                                                    
           do 260 i=km,k-1                                                      
             prd=0.d0                                                           
             do 270 j=1,idim                                                    
 270         prd=prd+v(j,i)*v(j,k)                                              
             do 280 j=1,idim                                                    
 280         v(j,k)=v(j,k)-prd*v(j,i)                                           
 260       continue                                                             
!                                                                               
           dnorm=0.0d0                                                          
           do 290 j=1,idim                                                      
 290       dnorm=dnorm+v(j,k)**2                                                
           s=1.d0/sqrt(dnorm)                                                   
           do 300 j=1,idim                                                      
 300       v(j,k)=v(j,k)*s                                                      
!                                                                               
        end if                                                                  
 250  continue                                                                  
      end                                                                       

!********* eigenvalues by the bisection method **********                       
!                                                                               
      subroutine bisec(alpha,beta,ndim,E,ne,eps)                                
!                                                                               
!    alpha  @ diagonal element                                                  
!    beta   @ subdiagonal element                                               
!    ndim   @ matrix dimension                                                  
!    E      # eigenvalues                                                       
!    ne     @ number of eigenvalues to calculate                                
!    eps    @ limit of error                                                    
                                                                                
      implicit real*8 (a-h,o-z)                                                 
      dimension alpha(ndim),beta(ndim),E(ne),b2(ndim)                           
!                                                                               
!      if(ndim.gt.2000)then                                                      
!          print *,' #(E04)# ndim given to bisec exceeds 2000'                   
!          stop                                                                  
!      end if                                                                    
!      if(ndim.gt.5000)then                                                      
!          print *,' #(E04)# ndim given to bisec exceeds 5000'                   
!          stop                                                                  
!      end if                                                                    

      if(ne.gt.ndim.or.ne.le.0)then                                             
          print *,' #(E05)# ne given to bisec out of range'                     
          stop                                                                  
      end if                                                                    
!                                                                               
!*** initial bound                                                              
      range=abs(alpha(1))+abs(beta(1))                                          
      do 10 k=2,ndim-1                                                          
 10   range=max(range,abs(beta(k-1))+abs(alpha(k))+abs(beta(k)))                
      range=max(range,abs(beta(ndim-1))+abs(alpha(ndim)))                       
      range=-range                                                              
!                                                                               
      b2(1)=0.d0                                                                
      do 20 i=2,ndim                                                            
 20   b2(i)=beta(i-1)**2                                                        
!                                                                               
      epsabs=abs(range)*eps                                                     
      do 30 i=1,ne                                                              
 30   E(i)=-range                                                               
      b=range                                                                   
!                                                                               
!*** bisection method                                                           
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
!                                                                               
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
