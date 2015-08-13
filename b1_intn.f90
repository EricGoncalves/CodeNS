module mod_b1_intn
  implicit none
contains
  subroutine b1_intn
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'intn'.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use kcle
    use schemanum
    use mod_convich
    use mod_mpi,only : rank
    implicit none
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=24) ::  cnumt
!
    call convich(knumt,cnumt)
!
    form='(/,2x,''initialisation du temps "numerique"'',/' &
         //'2x,''-----------------------------------'',/' &
         //'2x,''numt                     : '',11x,i5,2x,a)'
    if (rank==0) write(imp,form) numt,cnumt
!
    return
  end subroutine b1_intn
end module mod_b1_intn
