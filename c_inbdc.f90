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
    use boundary
    use sortiefichier
    use mod_b1_inbdc
    use mod_tcmd_inbdc
    use mod_inbdc
    implicit none
    integer          ::        iba, imot(nmx),       jba,       kba,     kibdc
    integer          ::        krr,     mfbea,     mfbeb, mnc(ip43),ncbd(ip41)
    integer          :: ncin(ip41),      nmot,l1,l2,mfbea1,mfbeb1
    double precision ::  epsmsh,   exs1,   exs2,x(ip21),y(ip21)
    double precision :: z(ip21)
!
!-----------------------------------------------------------------------
!
    character(len=2 ) :: tvi,tvj,tvk
    character(len=32) ::  mot(nmx)
!
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
