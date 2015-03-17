module mod_b1_svfw
  implicit none
contains
  subroutine b1_svfw(disc)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'svfw'.
!
!     INP
!_I    l      : arg int    ; numero de domaine
!_I    disc   : arg char   ; changement de discretisation (centre/noeud)
!_I    imp    : com int    ; unite logiq, sorties de controle
!
!***********************************************************************
!
    use sortiefichier
    implicit none
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=4 ) :: disc
!
    form='(/,2x,''sauvegarde des var_calc'',/' &
         //'2x,''-----------------------'',/' &
         //'2x,''discretisation           : '',12x,a)'
!
    write(imp,form) disc
!
    return
  end subroutine b1_svfw
end module mod_b1_svfw
