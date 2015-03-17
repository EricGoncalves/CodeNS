module mod_clsym
  implicit none
contains
  subroutine clsym( &
       mfb, &
       nxn,nyn,nzn,ncin,ncbd,v, &
       mmb,mpb,mpn)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement d'une condition de symetrie par rapport a la facette
!_A    frontiere.
!_A    Normales interieures.
!
!     INP
!_I    mfb        : arg int              ; numero de frontiere
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!
!     I/O
!_/    v          : arg real(ip11,ip60 ) ; variables a l'instant courant
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    implicit none
    integer          ::    m,  mb, mfb, mmb,  mn
    integer          ::  mpb, mpn,  mt,ncbd,ncin
    integer          ::   ni,  nl
    double precision ::  nxn, nyn, nzn, qn1, qtx
    double precision ::  qty, qtz, qx1, qy1, qz1
    double precision ::  ro1,roe1,   v
!
!-----------------------------------------------------------------------
!
    dimension v(ip11,ip60)
    dimension nxn(ip42),nyn(ip42),nzn(ip42)
    dimension ncin(ip41),ncbd(ip41)
    dimension mmb(mtt),mpb(mtt)
    dimension mpn(mtt)
!
    mt=mmb(mfb)
!
!!$OMP SIMD
    do m=1,mt
       mb  =mpb(mfb)+m
       mn  =mpn(mfb)+m
       nl  =ncbd(mb)
       ni  =ncin(mb)
!
       ro1 =v(ni,1)
       qx1 =v(ni,2)/ro1
       qy1 =v(ni,3)/ro1
       qz1 =v(ni,4)/ro1
       roe1=v(ni,5)
!
       qn1 =qx1*nxn(mn)+qy1*nyn(mn)+qz1*nzn(mn)
       qtx=qx1-qn1*nxn(mn)
       qty=qy1-qn1*nyn(mn)
       qtz=qz1-qn1*nzn(mn)
!
       v(nl,1)=ro1
       v(nl,2)=ro1*qtx
       v(nl,3)=ro1*qty
       v(nl,4)=ro1*qtz
       v(nl,5)=roe1
    enddo
!
    return
  end subroutine clsym
end module mod_clsym
