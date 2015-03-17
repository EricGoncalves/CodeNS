module mod_c_inbdc
  implicit none
contains
  subroutine c_inbdc( &
       mot,imot,nmot, &
       exs1,exs2, &
       x,y,z, &
       ncbd,ncin,mnc)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action inbdc.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use mod_b1_inbdc
    use mod_tcmd_inbdc
    use mod_inbdc
    implicit none
    integer          ::   iba, imot,  jba,  kba,kibdc
    integer          ::   krr,mfbea,mfbeb,  mnc, ncbd
    integer          ::  ncin, nmot
    double precision :: epsmsh,  exs1,  exs2,     x,     y
    double precision ::      z
!
!-----------------------------------------------------------------------
!
    character(len=2 ) :: tvi,tvj,tvk
    character(len=32) ::  mot(nmx)
!
    dimension imot(nmx)
    dimension x(ip21),y(ip21),z(ip21)
    dimension mnc(ip43)
    dimension ncbd(ip41)
    dimension ncin(ip41)
!
    call tcmd_inbdc( &
         mot,imot,nmot, &
         krr,mfbea,mfbeb,kibdc,epsmsh, &
         iba,jba,kba,tvi,tvj,tvk)
!
    if (kimp.ge.1) then
       call b1_inbdc( &
            krr,mfbea,mfbeb,kibdc,epsmsh, &
            iba,jba,kba,tvi,tvj,tvk)
    endif
!
    call inbdc( &
         exs1,exs2, &
         x,y,z, &
         ncbd,ncin,mnc, &
         krr,mfbea,mfbeb,kibdc,epsmsh, &
         iba,jba,kba,tvi,tvj,tvk)
!
    return
  end subroutine c_inbdc
end module mod_c_inbdc
