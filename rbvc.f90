      subroutine rbvc(t,ncbd,ncin,mnc)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables (t) sur des facettes frontieres coincidentes
!_A    (continuite par moyennes).
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
!     I/O
!_/    t          : arg real(ip11,ip60 ) ; variables de calcul
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use boundary
      use definition
implicit none
double precision :: t
integer :: ncbd
integer :: ncin
integer :: mnc
integer :: m
integer :: mb
integer :: mc
integer :: mf
integer :: mfb
integer :: mt
integer :: nc
integer :: nd
integer :: ndm
double precision :: tper
!
!-----------------------------------------------------------------------
!
      dimension t(ip11,ip60)
      dimension ncbd(ip41),ncin(ip41)
      dimension mnc(ip43)
!
      do mf=1,nbd
!
      mfb=lbd(mf)
      mt=mmb(mfb)
      tper=protat*float(mper(mfb))
!
      do m=1,mt
      mc=mpc(mfb)+m
      nc=mnc(mc)
      mb=mpb(mfb)+m
      nd=ncbd(mb)
      ndm=ncin(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
      t(nd,1) = 0.5*( t(ndm,1)+t(nc,1) )
      t(nd,2) = 0.5*( t(ndm,2)+t(nc,2) )
      t(nd,3) = 0.5*( t(ndm,3)+t(nc,3)*cos(tper)+t(nc,4)*sin(tper))
      t(nd,4) = 0.5*( t(ndm,4)+t(nc,4)*cos(tper)-t(nc,3)*sin(tper))
      t(nd,5) = 0.5*( t(ndm,5)+t(nc,5) )
!
      enddo
      enddo
!
      return
      end
