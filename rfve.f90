module mod_rfve
implicit none
contains
      subroutine rfve( &
                 t,ps,temp,cson, &
                 ncbd,ncin)
!
!***********************************************************************
!
!_DA    DATE fevrier 2007 - Eric Goncalves / LEGI
!
!     ACT
!_A    Calcul des variables (t) dans des mailles fictives adjacentes
!_A    a des frontieres, par extrapolation.
!
!_I    ip11       : arg int              ; dim, nbr max de cellules de tous les
!_I                                        dom (pts fictifs inclus)
!_I    ip41       : arg int              ; dim, nbr max de pts de ttes les front
!_I    ip60       : arg int              ; dim, nbr max d'equations
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use boundary
implicit none
double precision :: t
double precision :: ps
double precision :: temp
double precision :: cson
integer :: ncbd
integer :: ncin
integer :: m
integer :: mb
integer :: mf
integer :: mfb
integer :: mt
integer :: nd
integer :: ni
!
!-----------------------------------------------------------------------
!
      dimension t(ip11,ip60)
      dimension ncin(ip41),ncbd(ip41)
      dimension ps(ip11),temp(ip11),cson(ip11)
!
      do mf=1,nbd
       mfb=lbd(mf)
       mt=mmb(mfb)
!!$OMP SIMD
       do m=1,mt
        mb=mpb(mfb)+m
        nd=ncbd(mb)
        ni=ncin(mb)
        t(nd,1)=t(ni,1)
        t(nd,2)=t(ni,2)
        t(nd,3)=t(ni,3)
        t(nd,4)=t(ni,4)
        t(nd,5)=t(ni,5)
        ps(nd)=ps(ni)
        temp(nd)=temp(ni)
        cson(nd)=cson(ni)
       enddo
      enddo
!
      return
      end subroutine
end module
