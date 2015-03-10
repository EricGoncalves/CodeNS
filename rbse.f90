module mod_rbse
implicit none
contains
      subroutine rbse(t,ncin,ncbd)
!
!***********************************************************************
!
!     ACT
!_A    Calcul de la variable (t) sur des facettes frontieres par extrapolation.
!
!     INP
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
!     I/O
!_/    t          : arg real(ip11,ip60 ) ; variables de calcul
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
      dimension t(ip12)
      dimension ncbd(ip41),ncin(ip41)
!
!     definition de la variable aux bords (centre des facettes frontieres)
!
      do mf=1,nbd
!
      mfb=lbd(mf)
      mt=mmb(mfb)
!
!DEC$ IVDEP
      do m=1,mt
      ml=mpb(mfb)+m
      n=ncbd(ml)
      ni=ncin(ml)
      t(n) = t(ni)
      enddo
      enddo
!
      return
      end
end module
