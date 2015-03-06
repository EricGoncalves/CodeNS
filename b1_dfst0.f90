      subroutine b1_dfst0(roam,aam,tam)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfst0'.
!
!     INP
!_I    roam : arg real   ; etat de reference utilisateur
!_I                        dimensionne, masse volumique d'arret
!_I    aam  : arg real   ; etat de reference utilisateur
!_I                        dimensionne, vitesse du son d'arret
!_I    tam  : arg real   ; etat de reference utilisateur
!_I                        dimensionne, temperature d'arret
!_I    imp  : com int    ; unite logiq, sorties de controle
!
!***********************************************************************
!
      use sortiefichier
!
!-----------------------------------------------------------------------
!
      character *1316 form
!
           form='(/4x ,''grandeurs d''''arret amont :''/ &
             4x ,''~~~~~~~~~~~~~~~~~~~~~~~  ''/, &
             4x ,''roam   ='',e12.5,'' kg/m3'',/, &
             4x ,''aam    ='',e12.5,'' m/s  '',/, &
             4x ,''tam    ='',e12.5,'' K    '')'
      write(imp,form) roam,aam,tam
!
      return
      end
