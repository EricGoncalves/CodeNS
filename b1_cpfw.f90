      subroutine b1_cpfw
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'cpfw'.
!
!     INP
!_I    imp        : com int     ; unite logiq, sorties de controle
!
!***********************************************************************
!
      use sortiefichier
!
!-----------------------------------------------------------------------
!
      character *1316 form
!
       form='(/,2x,''realisation du calcul'',/' &
             //'2x,''---------------------'')'
      write(imp,form)
!
      return
      end
