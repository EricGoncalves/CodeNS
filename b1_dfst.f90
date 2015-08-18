module mod_b1_dfst
  implicit none
contains
  subroutine b1_dfst(nst)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dfst'.
!
!     INP
!_I    imp    : com int    ; unite logiq, sorties de controle
!_I    tnz    : com real   ; etat pour adimensionnement temperature
!_I    ronz   : com real   ; etat pour adimensionnement masse volumique
!_I    anz    : com real   ; etat pour adimensionnement vitesse du son d'arret
!_I    dnz    : com real   ; etat pour adimensionnement longueur
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_fige
    use sortiefichier
    use kcle
    use definition
    use mod_convich
    use mod_mpi,only:rank
    implicit none
    integer          :: nst
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
    character(len=24) ::  cd,cm,cal,cbe
    character(len=24) ::  ct,cro,ca
!
    call convich(kvarst(nst,4),cd)
    call convich(kvarst(nst,5),cm)
    call convich(kvarst(nst,6),cal)
    call convich(kvarst(nst,7),cbe)
    call convich(kvarst(nst,1),cro)
    call convich(kvarst(nst,2),ca)
    call convich(kvarst(nst,3),ct)
!
    if(rank==0) then
    form='(//,2x,''definition d''''un etat de type 1''/' &
         //'2x,''------------------------------'',/' &
         //'2x,''numero                   : '',9x,i2/' &
         //'2x,''masse volumique          : '',e11.4,''  kg/m3'',' &
         //'2x,a/' &
         //'2x,''vitesse du son           : '',e11.4,''    m/s'',' &
         //'2x,a/' &
         //'2x,''temperature              : '',e11.4,''      K'',' &
         //'2x,a/' &
         //'2x,''longueur                 : '',e11.4,''      m'',' &
         //'2x,a/' &
         //'2x,''nombre de Mach           : '',e11.4,''       '',' &
         //'2x,a/' &
         //'2x,''incidence                : '',e11.4,''    deg'',' &
         //'2x,a/' &
         //'2x,''derapage                 : '',e11.4,''    deg'',' &
         //'2x,a)'
    write(imp,form)nst, &
         varst(nst,1),cro, &
         varst(nst,2),ca, &
         varst(nst,3),ct, &
         varst(nst,4),cd, &
         varst(nst,5),cm, &
         varst(nst,6),cal, &
         varst(nst,7),cbe
    endif
!
    return
  end subroutine b1_dfst
end module mod_b1_dfst
