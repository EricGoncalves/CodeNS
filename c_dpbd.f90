module mod_c_dpbd
  implicit none
contains
  subroutine c_dpbd( &
       mot,imot,nmot, &
       ncbd,ncin, &
       nxn,nyn,nzn, &
       mnc, &
       mnr,xnr,ynr,znr, &
       tm1,tm2,tm3, &
       tm4,tm5,tm6)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dpbd.
!
!***********************************************************************
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use mod_tcmd_dpbd
    use mod_dpbd
    use mod_b1_dpbd
    implicit none
  integer          ::        img, imot(nmx),      lgrd,     lmfbd,        mf
  integer          ::       mfbe,        mg, mnc(ip43), mnr(ip44),ncbd(ip41)
  integer          :: ncin(ip41),      nmot
  double precision :: nxn(ip42),nyn(ip42),nzn(ip42),tm1(ip40),tm2(ip40)
  double precision :: tm3(ip40),tm4(ip40),tm5(ip40),tm6(ip40),xnr(ip44)
  double precision :: ynr(ip44),znr(ip44)
  integer         ,allocatable ::  lgr(:),lmfb(:)
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
    character(len=32) ::  typdat
 allocate(lmfb(mtb))
 allocate(lgr(nobj))
!
!
!
    call tcmd_dpbd( &
         mot,imot,nmot, &
         lmfb,lmfbd, &
         lgr,lgrd, &
         typdat)
!
    if (kimp.ge.1) then
       call b1_dpbd( &
            lmfb,lmfbd, &
            lgr,lgrd, &
            typdat)
    endif
!
    do mf=1,lmfbd
       mfbe=lmfb(mf)
       do mg=1,lgrd
          img=lgr(mg)
!
          call dpbd( &
               mfbe,img, &
               typdat, &
               ncbd,ncin, &
               nxn,nyn,nzn, &
               mnc, &
               mnr,xnr,ynr,znr, &
               tm1,tm2,tm3, &
               tm4,tm5,tm6)
!
       enddo
    enddo
!
    return
  end subroutine c_dpbd
end module mod_c_dpbd
