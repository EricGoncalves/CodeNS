module mod_clpara
implicit none
contains
      subroutine clpara( &
                 mfb, &
                 ncbd,v, &
                 mmb,mpb,ncin, &
                 pression,temp,cson)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement de la condition d'adherence
!_A    pour une paroi adiabatique.
!_A    Derivee de p : nulle suivant les lignes de maillage
!_A                   qui "sortent" de la paroi.
!
!
!_I    mfb        : arg int              ; numero de frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
implicit none
integer :: mfb
integer :: ncbd
double precision :: v
integer :: mmb
integer :: mpb
integer :: ncin
double precision :: pression
double precision :: temp
double precision :: cson
integer :: m
integer :: mb
integer :: mt
integer :: nc
integer :: ni
!
!-----------------------------------------------------------------------
!
      double precision rhoe1
      dimension v(ip11,ip60)
      dimension ncbd(ip41),ncin(ip41)
      dimension mmb(mtt),mpb(mtt)
      dimension pression(ip11),temp(ip11),cson(ip11)
!
      mt=mmb(mfb)
!
!!$OMP SIMD
      do m=1,mt
       mb=mpb(mfb)+m
       nc=ncbd(mb)
       ni=ncin(mb)
       rhoe1=v(ni,5)-0.5*(v(ni,2)**2+v(ni,3)**2+v(ni,4)**2)/v(ni,1)
!
       v(nc,1)=v(ni,1)
       v(nc,2)=0.
       v(nc,3)=0.
       v(nc,4)=0.
       v(nc,5)=rhoe1
!
       pression(nc)=pression(ni)      !hypothese dPdn=0
       temp(nc)=temp(ni)              !flux de chaleur nul: dTdn=0
       cson(nc)=cson(ni)
      enddo
!
      return
      end subroutine
end module
