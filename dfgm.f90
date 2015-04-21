module mod_dfgm
  implicit none
contains
  subroutine dfgm
!
!***********************************************************************
!
!     ACT
!_A    Calcul des valeurs de translation et de rotation periodiques pour
!_A    completer la definition du type de geometrie a calculer.
!
!     INP
!_I    config     : com char             ; type de config geometrique du calcul
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    pis2       : com real             ; pi divise par 2
!
!     OUT
!_O    ptrans     : com real             ; distance pour periodicite
!_O    protat     : com real             ; angle(rad) pour periodicite
!
!     I/O
!_/    perio      : com real             ; periodicite geometrique en angle ou
!_/                                        distance selon config
!
!***********************************************************************
!
    use constantes
    use definition
    use chainecarac
    implicit none
!
!-----------------------------------------------------------------------
!
!$OMP MASTER
    if ((config(1:3).eq.'gan').or.(config(1:3).eq.'hel')) then
       ptrans= 0.
       protat= 4.*pis2/perio
    else
       ptrans= perio
       protat= 0.
    endif
    if ((config(1:3).eq.'gan').or.(config(1:3).eq.'hel')) then
       perio=protat
    else
       perio=ptrans
    endif
!
!$OMP END MASTER
    return
  end subroutine dfgm
end module mod_dfgm
