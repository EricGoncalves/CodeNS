module mod_utprd
  implicit none
contains
  subroutine utprd( &
       bceqt, &
       mfl,pres,icyc, &
       mmb,mpb)
!
!***********************************************************************
!
!     ACT
!_A    Sous-programme utilisateur de preparation des donnees pour le
!_A    sous-programme clprd.
!_A    Il doit remplir le tableau pres des pressions.
!_A
!_I    mfl        : arg int              ; numero de la frontiere
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!
!     OUT
!_O    pres       : arg real(ip40      ) ; pression statique
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    implicit none
    integer          ::     icyc,       m,     mfl,      ml,mmb(mtt)
    integer          :: mpb(mtt),      mt
    double precision :: bceqt(ip41,neqt),      pres(ip40)
!
!-----------------------------------------------------------------------
!
!
    mt=mmb(mfl)
    do m=1,mt
       ml=mpb(mfl)+m
       pres(m)=bceqt(ml,1)
    enddo
!
  end subroutine utprd
end module mod_utprd
