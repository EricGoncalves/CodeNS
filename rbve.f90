      subroutine rbve(t,ncin,ncbd)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables (t) sur des facettes frontieres par extrapolation.
!
!     INP
!_I    ip11       : arg int              ; dim, nbr max de cellules de tous les
!_I                                        dom (pts fictifs inclus)
!_I    ip41       : arg int              ; dim, nbr max de pts de ttes les front
!_I    ip60       : arg int              ; dim, nbr max d'equations
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use boundary
implicit none
double precision :: t
integer :: ncin
integer :: ncbd
integer :: m
integer :: mf
integer :: mfb
integer :: ml
integer :: mt
integer :: n
integer :: ni
!
!-----------------------------------------------------------------------
!
      dimension t(ip11,ip60)
      dimension ncbd(ip41),ncin(ip41)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
      do mf=1,nbd
       mfb=lbd(mf)
       mt=mmb(mfb)
!DEC$ IVDEP
       do m=1,mt
        ml=mpb(mfb)+m
        n=ncbd(ml)
        ni=ncin(ml)
        t(n,1) = t(ni,1)
        t(n,2) = t(ni,2)
        t(n,3) = t(ni,3)
        t(n,4) = t(ni,4)
        t(n,5) = t(ni,5)
       enddo
      enddo
!
      return
      end
