module mod_rfsc
  implicit none
contains
  subroutine rfsc(t,ncbd,mnc)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables (t) dans des mailles fictives adjacentes a des
!_A    frontieres coincidentes (valeurs dans le domaine coincident).
!
!     INP
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
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
!_/    t          : arg real(ip11      ) ; variables de calcul
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    implicit none
    integer          ::    m,  mb,  mc,  mf, mfb
    integer          ::  mnc,  mt,  nc,ncbd,  nd
    double precision :: t
!
!-----------------------------------------------------------------------
!
    dimension t(ip11)
    dimension mnc(ip43),ncbd(ip41)
!
    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
!
!!$OMP SIMD
       do m=1,mt
          mc=mpc(mfb)+m
          nc=mnc(mc)
          mb=mpb(mfb)+m
          nd=ncbd(mb)
!
!     definition d'un scalaire aux points fictifs
!
          t(nd)=t(nc)
!
       enddo
    enddo
!
    return
  end subroutine rfsc
end module mod_rfsc
