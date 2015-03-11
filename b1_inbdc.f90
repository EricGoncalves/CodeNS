module mod_b1_inbdc
implicit none
contains
      subroutine b1_inbdc( &
                 krr,mfbea,mfbeb,kibdc,epsmsh, &
                 iba,jba,kba,tvi,tvj,tvk)
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'inbdc'.
!
!     INP
!_I    krr        : arg int              ; cle info sur front coinc
!_I    mfbea      : arg int              ; numero externe de frontiere
!_I    mfbeb      : arg int              ; numero externe de la front coinc
!_I    kibdc      : arg int              ; cle initialisation front coinc
!_I    epsmsh     : arg real             ; dist max entre deux pts confondus
!_I    iba        : arg int              ; ind i du pt coinc au pt d'indices
!_I                                        min de la frontiere
!_I    jba        : arg int              ; ind j du pt coinc au pt d'indices
!_I                                        min de la frontiere
!_I    kba        : arg int              ; ind k du pt coinc au pt d'indices
!_I                                        min de la frontiere
!_I    tvi        : arg char             ; sens de variation sur la front coinc
!_I                                        pour une var de l'ind i sur la front
!_I    tvj        : arg char             ; sens de variation sur la front coinc
!_I                                        pour une var de l'ind j sur la front
!_I    tvk        : arg char             ; sens de variation sur la front coinc
!_I                                        pour une var de l'ind k sur la front
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
      use para_fige
   use sortiefichier
implicit none
integer :: krr
integer :: mfbea
integer :: mfbeb
integer :: kibdc
double precision :: epsmsh
integer :: iba
integer :: jba
integer :: kba
!
!-----------------------------------------------------------------------
!
      character(len=1316) :: form
      character(len=2 ) :: tvi,tvj,tvk
!
          form='(/,2x,''calculs pour frontiere coincidente'',/' &
                //'2x,''----------------------------------'',/' &
                //'2x,''numero de la frontiere   : '',11x,i5/' &
                //'2x,''no frontiere qui coincide: '',11x,i5/' &
                //'2x,''cle calcul pts coincid.  : '',11x,i5)'
      write(imp,form) mfbea,mfbeb,kibdc
      if (kibdc.eq.1) then
          form='(2x,''cle positionnement front : '',11x,i5)'
      write(imp,form) krr
        if (krr.eq.1) then
          form='(2x,''dist max de 2 pts coinc. : '',4x,e12.5)'
      write(imp,form) epsmsh
        else if (krr.eq.0) then
          form='(2x,''indices du 1er pt coinc. : '',3(11x,i5)/' &
                //'2x,''correspondance des dir.  : '',3(14x,a))'
      write(imp,form) iba,jba,kba,tvi,tvj,tvk
        endif
      endif
!
      return
      end subroutine
end module
