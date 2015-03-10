      subroutine b1_dfgm
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfgm'.
!
!     INP
!_I    config     : com char             ; type de config geometrique du calcul
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    perio      : com real             ; periodicite geometrique en angle ou
!_I                                        distance selon config
!
!***********************************************************************
!
      use sortiefichier
   use definition
   use kcle
   use chainecarac
implicit none
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: form
      character(len=24) ::  cconfig,cperio
!
      call convich(kconfig,cconfig)
      call convich(kperio,cperio)
!
      if (config(1:3).eq.'gan') then
         form='(/,2x,''definition de la geometrie'',/' &
               //'2x,''--------------------------'',/' &
               //'2x,''config                   : '',12x,a,2x,a/' &
               //'2x,''nombre d''''aubes        : '',4x,e12.6,2x,a)'
      write(imp,form) config,cconfig, &
                      perio,cperio
      else if (config(1:3).eq.'hel') then
         form='(/,2x,''definition de la geometrie'',/' &
               //'2x,''--------------------------'',/' &
               //'2x,''config                   : '',12x,a,2x,a/' &
               //'2x,''nombre de pales          : '',4x,e12.6,2x,a)'
      write(imp,form) config,cconfig, &
                      perio,cperio
      else
         form='(/,2x,''definition de la geometrie'',/' &
               //'2x,''--------------------------'',/' &
               //'2x,''config                   : '',12x,a,2x,a/' &
               //'2x,''pas                      : '',4x,e12.6,2x,a)'
      write(imp,form) config,cconfig, &
                      perio,cperio
      endif
!
      return
      end
