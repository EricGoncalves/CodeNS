module mod_rbtc
implicit none
contains
      subroutine rbtc( &
                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                 ncbd,ncin,mnc)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des composantes du tenseur visqueux et des flux de chaleur
!_A    pour des facettes frontieres coincidentes (continuite par moyennes).
!
!     INP
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    mnc        : arg int (ip43      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule coincidente
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpc        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux front coinc
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!
!     OUT
!
!     I/O
!_/    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_/                                        contraintes visqueuses
!_/    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_/                                        contraintes visqueuses
!_/    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_/                                        contraintes visqueuses
!_/    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_/                                        contraintes visqueuses
!_/    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_/                                        contraintes visqueuses
!_/    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_/                                        contraintes visqueuses
!_/    qcx        : arg real(ip12      ) ; composante en x du flux de chaleur
!_/    qcy        : arg real(ip12      ) ; composante en y du flux de chaleur
!_/    qcz        : arg real(ip12      ) ; composante en z du flux de chaleur
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use boundary
      use definition
implicit none
double precision :: toxx
double precision :: toxy
double precision :: toxz
double precision :: toyy
double precision :: toyz
double precision :: tozz
double precision :: qcx
double precision :: qcy
double precision :: qcz
integer :: ncbd
integer :: ncin
integer :: mnc
double precision :: cr
integer :: m
integer :: mb
integer :: mc
integer :: mf
integer :: mfb
integer :: mt
integer :: nc
integer :: nd
integer :: ndm
double precision :: qcxr
double precision :: qcyr
double precision :: qczr
double precision :: sr
double precision :: txxr
double precision :: txyr
double precision :: txzr
double precision :: tyyr
double precision :: tyzr
double precision :: tzzr
!
!-----------------------------------------------------------------------
!
      dimension toxx(ip12),toxy(ip12),toxz(ip12), &
                toyy(ip12),toyz(ip12),tozz(ip12), &
                qcx(ip12),qcy(ip12),qcz(ip12)
      dimension ncin(ip41),ncbd(ip41)
      dimension mnc(ip43)
!
      do mf=1,nbd
!
      mfb=lbd(mf)
      mt=mmb(mfb)
      sr=-sin(real(mper(mfb))*protat)
      cr= cos(real(mper(mfb))*protat)
!
!!$OMP SIMD
      do m=1,mt
      mc =mpc(mfb)+m
      nc =mnc(mc)
      mb =mpb(mfb)+m
      nd =ncbd(mb)
      ndm=ncin(mb)
!
      txxr=toxx(nc)
      txyr=cr*toxy(nc)-sr*toxz(nc)
      txzr=cr*toxz(nc)+sr*toxy(nc)
      tyyr=cr*cr*toyy(nc)-2.*cr*sr*toyz(nc)+sr*sr*tozz(nc)
      tyzr=-cr*sr*(tozz(nc)-toyy(nc))+(2.*cr*cr-1.)*toyz(nc)
      tzzr=sr*sr*toyy(nc)+2.*cr*sr*toyz(nc)+cr*cr*tozz(nc)
      qcxr=qcx(nc)
      qcyr=cr*qcy(nc)-sr*qcz(nc)
      qczr=cr*qcz(nc)+sr*qcy(nc)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
      toxx(nd) = 0.5*( toxx(ndm)+txxr )
      toxy(nd) = 0.5*( toxy(ndm)+txyr )
      toxz(nd) = 0.5*( toxz(ndm)+txzr )
      toyy(nd) = 0.5*( toyy(ndm)+tyyr )
      toyz(nd) = 0.5*( toyz(ndm)+tyzr )
      tozz(nd) = 0.5*( tozz(ndm)+tzzr )
       qcx(nd) = 0.5*(  qcx(ndm)+ qcxr )
       qcy(nd) = 0.5*(  qcy(ndm)+ qcyr )
       qcz(nd) = 0.5*(  qcz(ndm)+ qczr )
!
      enddo
      enddo
!
      return
      end subroutine
end module
