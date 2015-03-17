module mod_clnrd
implicit none
contains
      subroutine clnrd( &
                 mfb, &
                 rod,roud,rovd,rowd,roed, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,l, &
                 pression,temp,cson)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement d'une condition de non reflexion.
!_A    Utilisation des cinq relations de compatibilite discretisees
!_A    combinees avec la donnee d'un etat de reference.
!_A      Dut=0
!_A      Dut=0
!_A      Dp - c2*Drho = 0
!_A      Dp + rho*c*Dun = 0
!_A      Dp - rho*c*Dun = 0
!_A    Normales interieures.
!
!     INP
!_I    mfb        : arg int              ; numero de frontiere
!_I    rod        : arg real(ip40      ) ; masse volumique donnee
!_I    roud       : arg real(ip40      ) ; ro*u donne, ro masse volumique
!_I                                        et u composante en x de la vitesse
!_I    rovd       : arg real(ip40      ) ; ro*v donne, ro masse volumique
!_I                                        et v composante en y de la vitesse
!_I    rowd       : arg real(ip40      ) ; ro*w donne, ro masse volumique
!_I                                        et w composante en z de la vitesse
!_I    roed       : arg real(ip40      ) ; ro*e donne, ro masse volumique
!_I                                        et e energie interne
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!
!     COM
!_C    Notation 0      : valeurs a l' instant n.
!_C    Notation s      : valeurs issues du schema.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use proprieteflu
implicit none
integer :: mfb
double precision :: rod
double precision :: roud
double precision :: rovd
double precision :: rowd
double precision :: roed
integer :: ncbd
double precision :: v
integer :: mmb
integer :: mpb
integer :: mpn
integer :: l
double precision :: pression
double precision :: temp
double precision :: cson
double precision :: am
double precision :: ap
double precision :: b0
double precision :: bs
double precision :: eps0
double precision :: epsm
double precision :: epsp
integer :: m
integer :: mb
integer :: mn
integer :: mt
integer :: n0c
integer :: n0n
integer :: nl
double precision :: p
double precision :: pd
double precision :: ps
double precision :: ps0
double precision :: qn
double precision :: qnd
double precision :: qns
double precision :: qtx
double precision :: qtxd
double precision :: qtxs
double precision :: qty
double precision :: qtyd
double precision :: qtys
double precision :: qtz
double precision :: qtzd
double precision :: qtzs
double precision :: qxd
double precision :: qxs
double precision :: qyd
double precision :: qys
double precision :: qzd
double precision :: qzs
double precision :: ro
double precision :: ro0
double precision :: roc0
double precision :: roqn0
double precision :: ros
!
!-----------------------------------------------------------------------
!
      double precision nxn,nyn,nzn
!
      dimension v(ip11,ip60)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
      dimension rod(ip40),roud(ip40),rovd(ip40),rowd(ip40),roed(ip40)
      dimension mmb(mtt),mpb(mtt)
      dimension mpn(mtt)
      dimension pression(ip11),temp(ip11),cson(ip11)
!
      n0n=npn(l)
      n0c=npc(l)
      mt=mmb(mfb)
!
!!$OMP SIMD
      do m=1,mt
       mb=mpb(mfb)+m
       mn=mpn(mfb)+m
       nl=ncbd(mb)
!      etat '0' au pas de temps n
       ro0=v(nl,1)
       ps0=gam1*(v(nl,5)-0.5*(v(nl,2)**2+v(nl,3)**2+v(nl,4)**2)/v(nl,1)-pinfl)
       roc0=sqrt(gam*ps0*ro0)
       roqn0=v(nl,2)*nxn(mn)+v(nl,3)*nyn(mn)+v(nl,4)*nzn(mn)
       epsm=.5+sign(.5, roc0-roqn0)
       eps0=.5+sign(.5,-roqn0)
       epsp=.5+sign(.5,-roc0-roqn0)
!      etat de reference (valeur prise a t=0)
       qxd=roud(m)/rod(m)
       qyd=rovd(m)/rod(m)
       qzd=rowd(m)/rod(m)
       pd=gam1*(roed(m)-.5*rod(m)*(qxd**2+qyd**2+qzd**2)-pinfl)
       qnd=qxd*nxn(mn)+qyd*nyn(mn)+qzd*nzn(mn)
       qtxd=qxd-qnd*nxn(mn)
       qtyd=qyd-qnd*nyn(mn)
       qtzd=qzd-qnd*nzn(mn)
!      etat schema
       ros=v(nl,1)
       qxs=v(nl,2)/ros
       qys=v(nl,3)/ros
       qzs=v(nl,4)/ros
       ps=gam1*(v(nl,5)-.5*ros*(qxs**2+qys**2+qzs**2)-pinfl)
       qns=qxs*nxn(mn)+qys*nyn(mn)+qzs*nzn(mn)
       qtxs=qxs-qns*nxn(mn)
       qtys=qys-qns*nyn(mn)
       qtzs=qzs-qns*nzn(mn)
!
       qtx=eps0*qtxs         +(1.-eps0)*qtxd
       qty=eps0*qtys         +(1.-eps0)*qtyd
       qtz=eps0*qtzs         +(1.-eps0)*qtzd
       am=epsm*(ps-roc0*qns)+(1.-epsm)*(pd-roc0*qnd)
       ap=epsp*(ps+roc0*qns)+(1.-epsp)*(pd+roc0*qnd)
       qn=0.5*(ap-am)/roc0
       p=0.5*(ap+am)
       bs=(p-ps)*ro0**2/roc0**2+ros
       b0=(p-pd)*ro0**2/roc0**2+rod(m)
       ro=eps0*bs+(1.-eps0)*b0
!
       v(nl,1)=ro
       v(nl,2)=ro*(qtx+qn*nxn(mn))
       v(nl,3)=ro*(qty+qn*nyn(mn))
       v(nl,4)=ro*(qtz+qn*nzn(mn))
       v(nl,5)=p/gam1+pinfl+0.5*(v(nl,2)**2+v(nl,3)**2+v(nl,4)**2)/ro
!
       pression(nl)=p
       temp(nl)=gam*pression(nl)/ro
       cson(nl)=sqrt(temp(nl))
      enddo
!
      return
      end subroutine
end module
