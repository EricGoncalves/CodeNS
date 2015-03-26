module mod_c_infw
  implicit none
contains
  subroutine c_infw( &
       mot,imot,nmot, &
       x,y,z,v,mut,tnte1,utau, &
       vdual,vdual1,vdual2)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action infw.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use modeleturb
    use mod_infw
    use mod_tcmd_infw
    use mod_b1_infw
    implicit none
  integer          :: imot(nmx),     kina,        l,    ldomd,     nmot
  double precision ::         mut(ip12), tnte1(ip11,ip60),       utau(ip42),     v(ip11,ip60), vdual(ip11,ip60)
  double precision :: vdual1(ip11,ip60),vdual2(ip11,ip60),          x(ip21),          y(ip21),          z(ip21)
  integer         ,allocatable :: ldom(:)
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
 allocate(ldom(nobj))
!
!
    call tcmd_infw( &
         mot,imot,nmot, &
         ldom,ldomd,kina)
!
    if (kimp.ge.1) then
       call b1_infw(ldom,ldomd,kina)
    endif
!
    do l=1,ldomd
!
       call infw( &
            ldom(l),x,y,z,v,mut,tnte1, &
            kina,utau, &
            vdual,vdual1,vdual2)
!
    enddo
!
    return
  end subroutine c_infw
end module mod_c_infw
