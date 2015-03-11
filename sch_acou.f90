module mod_sch_acou
implicit none
contains
      subroutine sch_acou( &
                 lm,ityprk, &
                 u,ff,dtpas, &
                 x,y,z,idcyc, &
                 vol)
!
!***********************************************************************
!
!_DA  DATE_C : avril 2002 - Eric Goncalves / SINUMEF
!
!     ACT
!_A    Bilan de flux avec source acoustique.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use proprieteflu
      use maillage
      use constantes
      use schemanum
implicit none
integer :: indc
integer :: i
integer :: j
integer :: k
integer :: lm
integer :: ityprk
double precision :: u
double precision :: ff
double precision :: dtpas
double precision :: x
double precision :: y
double precision :: z
integer :: idcyc
double precision :: vol
integer :: i1
integer :: i2
integer :: i2m1
integer :: ind1
integer :: ind2
integer :: j1
integer :: j2
integer :: j2m1
integer :: k1
integer :: k2
integer :: k2m1
integer :: m
integer :: n0c
integer :: nc
integer :: nid
integer :: nijd
integer :: njd
double precision :: omeg
double precision :: temps
double precision :: ts1
double precision :: xcc
double precision :: ycc
double precision :: zcc
!
!-----------------------------------------------------------------------
!
      dimension u(ip11,ip60),ff(ip11,ip60)
      dimension vol(ip11)
      dimension x(ip21),y(ip21),z(ip21)
!
      indc(i,j,k)=1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
!
      n0c=npc(lm)
      i1=ii1(lm)
      i2=ii2(lm)
      j1=jj1(lm)
      j2=jj2(lm)
      k1=kk1(lm)
      k2=kk2(lm)
!
      nid = id2(lm)-id1(lm)+1
      njd = jd2(lm)-jd1(lm)+1
      nijd = nid*njd
!
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
!-----calcul des accroissements par maille :
!
      ind1 = indc(i1  ,j1  ,k1  )
      ind2 = indc(i2m1,j2m1,k2m1)
!     pulsation
      omeg=4.*pis2*freq
!
      if(ityprk.eq.0) then
       do m=ind1,ind2
        nc=m+n0c
!       coordonnees centre facette
        xcc=(x(nc)     +x(nc     +1)+x(nc     +nid)+x(nc     +nid+1) &
            +x(nc+nijd)+x(nc+nijd+1)+x(nc+nijd+nid)+x(nc+nijd+nid+1))* &
             0.125
        ycc=(y(nc)     +y(nc     +1)+y(nc     +nid)+y(nc     +nid+1) &
            +y(nc+nijd)+y(nc+nijd+1)+y(nc+nijd+nid)+y(nc+nijd+nid+1))* &
             0.125
        zcc=(z(nc)     +z(nc     +1)+z(nc     +nid)+z(nc     +nid+1) &
            +z(nc+nijd)+z(nc+nijd+1)+z(nc+nijd+nid)+z(nc+nijd+nid+1))* &
             0.125
!
        temps=real(idcyc)*dtpas
!        ts1=0.5*exp(-log(2.)*((xcc-x0)**2 + (ycc-y0)**2
!     &           +  (zcc-z0)**2)/cga**2)*cos(omeg*temps)
        ts1=0.5*exp(-log(2.)*((xcc-x0)**2+(ycc-y0)**2)/cga**2) &
                 *cos(omeg*temps)
        u(nc,1) = u(nc,1) - ts1*vol(nc)
        u(nc,5) = u(nc,5) - ts1*vol(nc)/(gam-1.)
       enddo
!
      else
!
       do m=ind1,ind2
        nc=m+n0c
!       coordonnees centre facette
        xcc=(x(nc)     +x(nc     +1)+x(nc     +nid)+x(nc     +nid+1) &
            +x(nc+nijd)+x(nc+nijd+1)+x(nc+nijd+nid)+x(nc+nijd+nid+1))* &
             0.125
        ycc=(y(nc)     +y(nc     +1)+y(nc     +nid)+y(nc     +nid+1) &
            +y(nc+nijd)+y(nc+nijd+1)+y(nc+nijd+nid)+y(nc+nijd+nid+1))* &
             0.125
        zcc=(z(nc)     +z(nc     +1)+z(nc     +nid)+z(nc     +nid+1) &
            +z(nc+nijd)+z(nc+nijd+1)+z(nc+nijd+nid)+z(nc+nijd+nid+1))* &
             0.125
!
        temps=real(idcyc)*dtpas
!        ts1=1.*exp(-log(2.)*((xcc-x0)**2 + (ycc-y0)**2
!     &           +  (zcc-z0)**2)/cga**2)*cos(omeg*temps)
        ts1=0.5*exp(-log(2.)*((xcc-x0)**2+(ycc-y0)**2)/cga**2) &
                 *cos(omeg*temps)
        u(nc,1) = u(nc,1) + ff(nc,1) - ts1*vol(nc)
        u(nc,5) = u(nc,5) + ff(nc,5) - ts1*vol(nc)/(gam-1.)
       enddo
      endif
!
      return
      end subroutine
end module
