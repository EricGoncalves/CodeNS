module mod_utpari
  implicit none
contains
  subroutine utpari( &
       bceqt, &
       mfl,tp, &
       mmb,mpb)
!
!***********************************************************************
!
!     ACT
!_A    Sous-programme utilisateur de preparation des donnees pour le
!_A    sous-programme clpari.
!_A    Il doit remplir le tableau tp de la temperature paroi.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    implicit none
    integer          ::        m,     mfl,      ml,mmb(mtt),mpb(mtt)
    integer          ::       mt
    double precision :: bceqt(ip41,neqt),        tp(ip40)
!
!-----------------------------------------------------------------------
!
!
    mt=mmb(mfl)
!$OMP SIMD
    do m=1,mt
       ml=mpb(mfl)+m
       tp(m)=bceqt(ml,1)
    enddo
!
    return
  end subroutine utpari
end module mod_utpari
