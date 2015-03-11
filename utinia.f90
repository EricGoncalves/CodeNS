module mod_utinia
implicit none
contains
      subroutine utinia( &
                 l,x,y,z,v,mut, &
                 kina, &
                 vdual,vdual1,vdual2)
!
!***********************************************************************
!
!     ACT
!         Lecture du fichier finia1
!         Initialisation des variables conservatives
!         La viscosite turbulente mut est mise a 0
!         Initialisation de vdual pour pas de temps dual
!
!     INP
!_I    out        : com int              ; unite logiq, moyennes des residus
!_I    sec        : com int              ; unite logiq, pts a p ou ro negatifs
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use sortiefichier
      use chainecarac
      use maillage  
      use definition
      use proprieteflu
      use schemanum
implicit none
integer :: inc
integer :: ind
integer :: id
integer :: jd
integer :: kd
integer :: i
integer :: j
integer :: k
integer :: l
double precision :: x
double precision :: y
double precision :: z
double precision :: v
integer :: kina
double precision :: vdual
double precision :: vdual1
double precision :: vdual2
double precision :: a
double precision :: alpha
double precision :: alphar
double precision :: beta
double precision :: betar
double precision :: degrad
integer :: i1
integer :: i2
integer :: i2m1
integer :: j1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k2
integer :: k2m1
integer :: n
integer :: n0
integer :: nci
integer :: ncij
integer :: ncijk
integer :: ncik
integer :: ncj
integer :: ncjk
integer :: nck
integer :: nid
integer :: nijd
integer :: njd
double precision :: p
double precision :: pis2
double precision :: q
double precision :: rmach
double precision :: ro
double precision :: rou
double precision :: rov
double precision :: row
double precision :: ym
double precision :: zm
!
!-----------------------------------------------------------------------
!
      real mut
      dimension x(*),y(*),z(*)
      dimension v(ip11,ip60)
      dimension mut(*)
      dimension vdual(ip11,ip60),vdual1(ip11,ip60),vdual2(ip11,ip60)
!
      inc(id,jd,kd)=id+jd*nid+kd*nijd
      ind(i,j,k)=n0+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
!
      n0=npc(l)
      i1=ii1(l)
      i2=ii2(l)
      j1=jj1(l)
      j2=jj2(l)
      k1=kk1(l)
      k2=kk2(l)
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd= nid*njd
!
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
      nci  =inc(1,0,0)
      ncj  =inc(0,1,0)
      nck  =inc(0,0,1)
      ncij =inc(1,1,0)
      ncik =inc(1,0,1)
      ncjk =inc(0,1,1)
      ncijk=inc(1,1,1)
!
      pis2=atan2(1.,0.)
      degrad=pis2/90.
!
!----------if en commentaire pour cas multidomaines
!     if(l.eq.1) then
        write(imp,'(/,"===>utinia: ouverture et lecture inia1")')
        close(inia1)
        open(inia1,file='finia1')
!
        read(inia1,*) rmach,alpha,beta
        close(inia1)
        write(imp,'(12x,"rmach=",f10.3,4x,"alpha=",f10.3,4x,"beta=",f10.4)')rmach,alpha,beta
        vrtmac=rmach
        vrtalp=alpha
        if(abs(beta).le.tiny(1.)) then
          write(imp,'(/,"!!!utinia: probleme possible avec la ",' &
          //'"correction de vorticite. alpha doit etre l incidence.",/,' &
          //'10x,"le profil doit etre en x et z")')
        end if
!     endif
!
!      print*,'===>utinia:  mach=',rmach,' alpha=',alpha,' beta=',beta
!
      ro=roa1/(1.+gam2*rmach**2)**gam4
      a=aa1/(1.+gam2*rmach**2)**.5
!      p=pa1/(1.+gam2*rmach**2)**(gam/gam1)
      p=(1.+gam2*rmach**2)**(-gam/gam1)/gam
      q=rmach*a
      alphar=alpha*degrad
      betar=beta*degrad
      rou=ro*q*cos(alphar)*cos(betar)
      rov=-ro*q*sin(betar)
      row=ro*q*sin(alphar)*cos(betar)
!
      do k=k1,k2m1
       do j=j1,j2m1
        do i=i1,i2m1
         n=ind(i,j,k)
         ym = 0.125*( y (n     )+y (n+nci  ) &
                     +y (n+ncj )+y (n+ncij ) &
                     +y (n+nck )+y (n+ncik ) &
                     +y (n+ncjk)+y (n+ncijk) )
         zm = 0.125*( z (n     )+z (n+nci  ) &
                     +z (n+ncj )+z (n+ncij ) &
                     +z (n+nck )+z (n+ncik ) &
                     +z (n+ncjk)+z (n+ncijk) )
!
         v(n,1)=ro
         v(n,2)=rou
         v(n,3)=rov + (ro*zm*omg)
         v(n,4)=row - (ro*ym*omg)
         v(n,5)=p/gam1 + pinfl + 0.5*(v(n,2)**2+v(n,3)**2+v(n,4)**2)/ro
        enddo
       enddo
      enddo
!
      if(equat(1:2).eq.'ns') then
       do k=k1,k2m1
        do j=j1,j2m1
         do i=i1,i2m1
          n=ind(i,j,k)
          mut(n)=0.
         enddo
        enddo
       enddo
      endif
!
      if(kfmg.eq.3) then
       do k=k1,k2m1
        do j=j1,j2m1
          do i=i1,i2m1
           n = ind(i,j,k)
           vdual(n,1) = v(n,1)
           vdual(n,2) = v(n,2)
           vdual(n,3) = v(n,3)
           vdual(n,4) = v(n,4)
           vdual(n,5) = v(n,5)
!
           vdual1(n,1) = v(n,1)
           vdual1(n,2) = v(n,2)
           vdual1(n,3) = v(n,3)
           vdual1(n,4) = v(n,4)
           vdual1(n,5) = v(n,5)
!
           vdual2(n,1) = v(n,1)
           vdual2(n,2) = v(n,2)
           vdual2(n,3) = v(n,3)
           vdual2(n,4) = v(n,4)
           vdual2(n,5) = v(n,5)
          enddo
         enddo
        enddo
      endif
!
      return
      end subroutine
end module
