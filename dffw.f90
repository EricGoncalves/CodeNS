      subroutine dffw
!
!***********************************************************************
!
!     ACT
!_A    Calcul de la vitesse de rotation du repere relatif en rad.s-1 pour
!_A    completer la definition du type d'ecoulement a calculer.
!
!     INP
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    klomg      : com int              ; cle pour rotation du repere relatif
!_I    pis2       : com real             ; pi divise par 2
!
!     I/O
!_/    omg        : com real             ; vitesse rotation du repere relatif
!
!***********************************************************************
!
      use constantes
      use definition
      use chainecarac
   use maillage
implicit none
!
!-----------------------------------------------------------------------
!
      omg=4.*pis2/60.* omg
      omg=omg*dnz/anz
!
      neqtx=5
      if (equat(6:7).eq.'ke') then
       neqtx=7
      endif
!
      return
      end
