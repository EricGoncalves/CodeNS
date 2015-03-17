module mod_smg_cf
  implicit none
contains
  subroutine smg_cf( &
       imgc,imgf, &
       vol , &
       vv,vc)
!***********************************************************************
!
!     ACT
!_A    Coarse to Fine : Correction from Coarse --> Fine
!_A    Prolongation
!_A                       COARSE level : imgc
!_A                       FINE   level : imgf
!
!-----------------------------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use chainecarac
    use schemanum
    use mod_smg_cf_2d
    use mod_smg_cf_3d
    implicit none
    integer          :: imgc,imgf
    double precision ::  vc,vol, vv
!
!-----------------------------------------------------------------------
!
    dimension vol (ip11)
    dimension vv(ip11,ip60),vc(ip11,ip60)
!
    ktrans=1
    if(equat(3:5).eq.'2dk' .or. equat(3:5).eq.'2xk' ) then
!                       ----------      --------------------    ---------
! ---    Interpolation  DISTANCES AUX NOEUDS ou BILINEAIRE   ou VOLUMIQUE  ---
!                       ----------      --------------------    ---------
       if ( (ktrans.eq.1).or.(ktrans.eq.2).or.(ktrans.eq.3) ) then
          call smg_cf_2d ( &
               imgc,imgf, &
               vol, &
               vv,vc)
       else
          write(6,*) ' ********************************************* '
          write(6,*) ' Cle Operateurs de transferts Coarse-Fine '
          write(6,*) ' ********************************************* '
          stop
       endif
!
    elseif(equat(3:5).eq.'2dj') then
       write(6,*) ' ********************************************* '
       write(6,*) '          Pas de multigrille en 2dj            '
       write(6,*) ' ********************************************* '
       stop
!
    elseif(equat(3:5).eq.'2di') then
       write(6,*) ' ********************************************* '
       write(6,*) '          Pas de multigrille en 2di            '
       write(6,*) ' ********************************************* '
       stop
!
    elseif(equat(3:4).eq.'3d') then
       call smg_cf_3d ( &
            imgc,imgf, &
            vol, &
            vv,vc)
!
    else
       stop 'smg_cf'
    end if
!
    return
  end subroutine smg_cf
end module mod_smg_cf
