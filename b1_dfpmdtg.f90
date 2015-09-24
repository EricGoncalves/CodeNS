module mod_b1_dfpmdtg
  implicit none
contains
  subroutine b1_dfpmdtg
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfpmdtdg'.
!
!     INP
!_I    kdtl       : arg int    ; cle d'utilisation pas de temps local
!_I    icychr0    : arg int  ; nbr de cycl en deb de calcul au cours
!_I                            desquelles le pas de temps mis a jour
!_I                            a chaque it
!_I    ncychro    : arg int   ; freq en it de mise a jour dupas de temps
!_I    dt1min     : arg real  ; pas de temps constant desire
!_I    imp        : com int   ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use schemanum
    use kcle
    use mod_convich
    implicit none
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=24) ::  ckdtl,cicychr0,cncychro,cdt1min
!
    call convich(kkdtl,ckdtl)
    call convich(kicychr0,cicychr0)
    call convich(kncychro,cncychro)
    call convich(kdt1min,cdt1min)
!
    form='(/,2x,''avance en temps'',/' &
         //'2x,''---------------'',/' &
         //'2x,''kdtl                     : '',11x,i5,2x,a/' &
         //'2x,''icychr0                  : '',11x,i5,2x,a/' &
         //'2x,''ncychro                  : '',11x,i5,2x,a/' &
         //'2x,''dt1min                   : '',4x,e12.6,2x,a)'
!
    write(imp,form) kdtl,ckdtl, &
         icychr0,cicychr0, &
         ncychro,cncychro, &
         dt1min,cdt1min
!
    return
  end subroutine b1_dfpmdtg
end module mod_b1_dfpmdtg
