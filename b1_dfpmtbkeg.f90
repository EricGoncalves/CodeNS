      subroutine b1_dfpmtbkeg
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfpmtbkeg'.
!
!     INP
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    gam        : com real             ; rapport des chaleurs specifiques
!_I    pr         : com real             ; nombre de Prandtl
!_I    prt        : com real             ; nombre de Prandtl turbulent
!_I    reynz      : com real             ; nombre de Reynolds calcule avec
!_I                                        les grandeurs d'adimensionnement,
!_I                                        pour definir la loi de Sutherland
!
!-----------------------------------------------------------------------
!
      use sortiefichier
      use modeleturb
      use schemanum
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: form
!
       form='(/,2x,''definition des parametres ke''/' &
             //'2x,''----------------------------'',/' &
             //'2x,''rokinf                   : '',4x,e12.6/' &
             //'2x,''roeinf                   : '',4x,e12.6/' &
             //'2x,''epsk                     : '',4x,e12.6/' &
             //'2x,''epse                     : '',4x,e12.6/' &
             //'2x,''cke1                     : '',4x,e12.6/' &
             //'2x,''cke2                     : '',4x,e12.6/' &
             //'2x,''alfak                    : '',4x,e12.6/' &
             //'2x,''alfae                    : '',4x,e12.6/' &
             //'2x,''rki2t                    : '',4x,e12.6/' &
             //'2x,''rki4t                    : '',4x,e12.6)'
      write(imp,form) rokinf, &
                      roeinf, &
                      epsk, &
                      epse, &
                      cke1, &
                      cke2, &
                      alfak, &
                      alfae, &
                      rki2t, &
                      rki4t  
!
      return
      end
