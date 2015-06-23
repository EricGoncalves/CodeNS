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
    do l1=1,mtb
      if (new2old_f(l1)==mfbea) then
        mfbea1=l1
        do l2=1,mtb
          if (new2old_f(l1)==mfbeb) then
            mfbeb1=l2

!                      imin=min(max(ii1(l2),save_iminb(fr)),ii2(l2))
!                      imax=min(max(ii1(l2),save_imaxb(fr)),ii2(l2))
!                      jmin=min(max(jj1(l2),save_jminb(fr)),jj2(l2))
!                      jmax=min(max(jj1(l2),save_jmaxb(fr)),jj2(l2))
!                      kmin=min(max(kk1(l2),save_kminb(fr)),kk2(l2))
!                      kmax=min(max(kk1(l2),save_kmaxb(fr)),kk2(l2))
!                      if ((imax-imin>=1.and.jmax-jmin>=1).or. &
!                           (imax-imin>=1.and.kmax-kmin>=1).or. &
!                           (jmax-jmin>=1.and.kmax-kmin>=1)) then

            call inbdc( &
                 exs1,exs2, &
                 x,y,z, &
                 ncbd,ncin,mnc, &
                 krr,mfbea,mfbeb,kibdc,epsmsh, &
                 iba,jba,kba,tvi,tvj,tvk)
          endif
        enddo
      endif
    enddo
!
    return
  end subroutine c_inbdc
end module mod_c_inbdc
