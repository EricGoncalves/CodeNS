      subroutine b1_end
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'end'.
!
!***********************************************************************
!
      use sortiefichier
!
!-----------------------------------------------------------------------
!
      character *1316 form
!
          form='(/,2x,''essai termine'',/' &
                //'2x,''-------------'',)'
      write(imp,form)
!
      return
      end
