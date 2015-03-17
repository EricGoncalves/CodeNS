module mod_b1_cpbd
  implicit none
contains
  subroutine b1_cpbd
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'cpbd'.
!
!     INP
!_I    kexl       : arg int              ; cle traitement frontieres
!_I                                        avant exploitation
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!
    use sortiefichier
    use boundary
    use kcle
    use mod_convich
    implicit none
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=24) :: ckexl
!
    call convich(kkexl,ckexl)
!
    form='(/,2x,''application de conditions aux limites'',/' &
         //'2x,''-------------------------------------''/' &
         //'2x,''cle d''''initialisation     : '',11x,i5,2x,a)'
    write(imp,form) kexl,ckexl
!
    return
  end subroutine b1_cpbd
end module mod_b1_cpbd
