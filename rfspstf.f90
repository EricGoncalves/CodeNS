module mod_rfspstf
  implicit none
contains
  subroutine rfspstf(t0,ncin,ncbd)
!
!***********************************************************************
!
!     ACT
!_A    Determination des variables de calcul aux centres de mailles fictives
!_A    par extrapolation a partir des valeurs aux centres des mailles
!_A    interieures jouxtant les frontieres et des valeurs sur les facettes
!_A    frontieres.
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
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    implicit none
    integer          ::    m,  mf, mfb,  ml,  mt
    integer          ::    n,ncbd,ncin
    double precision :: t0
!
!-----------------------------------------------------------------------
!
    dimension t0(ip11)
    dimension ncin(ip41),ncbd(ip41)
!
!     definition des variables aux points fictifs
!
    do mf=1,nbd
!
       mfb=lbd(mf)
       mt=mmb(mfb)
!
!!!$OMP SIMD
       do m=1,mt
          ml=mpb(mfb)+m
          n=ncbd(ml)
          t0(n) = 2*t0(n)-t0(ncin(ml))
       enddo
    enddo
!
    return
  end subroutine rfspstf
end module mod_rfspstf
