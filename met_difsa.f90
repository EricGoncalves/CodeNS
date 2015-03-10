module mod_met_difsa
implicit none
contains
      subroutine met_difsa( &
                 l, &
                 sn, &
                 vol,s, &
                 t,dtdx,dtdy,dtdz, &
                 txxf5x,txyf5y,txzf5z,sdif, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!      Calcul de grad(nutilde).grad(rho*nutilde)
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    s          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!
!     OUT
!-O    sdif       : arg real(ip00    )  ; cb2/sigma (grad(nu tilde).grad(
!_O                                       ro nu tilde))
!
!     I/O
!_/    txxf5x     : arg real(ip12     )  ; comp x grad(nu_tilde) puis
!_/                                        tenseur visqueux
!_/    txyf5y     : arg real(ip12     )  ; comp y grad(nu_tilde) puis
!_/                                        tenseur visqueux
!_/    txzf5z     : arg real(ip12     )  ; comp z grad(nu_tilde) puis
!_/                                        tenseur visqueux
!
!     LOC
!_L    t          : arg real(ip00     )  ; variable de travail
!_L    dtdx       : arg real(ip00     )  ; grad(t)  t,x
!_L    dtdy       : arg real(ip00     )  ; grad(t)  t,y
!_L    dtdz       : arg real(ip00     )  ; grad(t)  t,z
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use modeleturb
      use chainecarac
use mod_teq_grads
implicit none
integer :: indc
integer :: i
integer :: j
integer :: k
integer :: l
double precision :: sn
double precision :: vol
double precision :: s
double precision :: t
double precision :: dtdx
double precision :: dtdy
double precision :: dtdz
double precision :: txxf5x
double precision :: txyf5y
double precision :: txzf5z
double precision :: sdif
double precision :: cmui1
double precision :: cmui2
double precision :: cmuj1
double precision :: cmuj2
double precision :: cmuk1
double precision :: cmuk2
integer :: i1
integer :: i1m1
integer :: i2
integer :: i2m1
integer :: imax
integer :: imin
integer :: ind1
integer :: ind2
integer :: j1
integer :: j1m1
integer :: j2
integer :: j2m1
integer :: jmax
integer :: jmin
integer :: k1
integer :: k1m1
integer :: k2
integer :: k2m1
integer :: kmax
integer :: kmin
integer :: lgsnlt
integer :: m
integer :: n
integer :: n0c
integer :: nid
integer :: nijd
integer :: njd
integer :: npsn
!
!-----------------------------------------------------------------------
!
      dimension s(ip11,ip60)
      dimension sn(ip31*ndir), &
                vol(ip11)
      dimension t(ip00),dtdx(ip00),dtdy(ip00),dtdz(ip00),sdif(ip00)
      dimension txxf5x(ip12),txyf5y(ip12),txzf5z(ip12)
      dimension cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
      indc(i,j,k)=n0c+1+(i-id1(l))+(j-jd1(l))*nid+(k-kd1(l))*nijd
!
      n0c=npc(l)
      i1=ii1(l)
      i2=ii2(l)
      j1=jj1(l)
      j2=jj2(l)
      k1=kk1(l)
      k2=kk2(l)
!
      i1m1=i1-1
      j1m1=j1-1
      k1m1=k1-1
!
      i2m1=i2-1
      j2m1=j2-1
      k2m1=k2-1
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd= nid*njd
!
      imin=i1m1
      imax=i2
      jmin=j1m1
      jmax=j2
      kmin=k1m1
      kmax=k2
!
!calcul de grad(ro nu-tilde))
!
      if (equat(3:5).eq.'2di') then
         imin=i1
         imax=i2m1
      endif
      if (equat(3:5).eq.'2dj') then
         jmin=j1
         jmax=j2m1
      endif
      if (equat(3:5).eq.'2dk') then
         kmin=k1
         kmax=k2m1
      endif
!
      do k=kmin,kmax
         do j=jmin,jmax
            ind1=indc(imin,j,k)
            ind2=indc(imax,j,k)
            do n=ind1,ind2
               m=n-n0c
               t (m)= s(n,6)
            enddo
         enddo
      enddo
!
      npsn  =ndir*npfb(l)+1
      lgsnlt=nnn(l)
!
      call teq_grads( &
                 l, &
                 equat, &
                 sn(npsn),lgsnlt, &
                 vol, &
                 t   , &
                 dtdx,dtdy,dtdz, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)   !cEG-TOM
!
! calcul de cb2/sigma (grad(nu tilde).grad(ro nu tilde))
!
      ind1=indc(i1  ,j1  ,k1  )
      ind2=indc(i2m1,j2m1,k2m1)
      do n=ind1,ind2
         m=n-n0c
         sdif(m)=dtdx(m)*txxf5x(n)+dtdy(m)*txyf5y(n)+ &
                 dtdz(m)*txzf5z(n)
         sdif(m)=sdif(m)*cb2/sigma
      enddo
      return
      end
end module
