module mod_dpbd
  implicit none
contains
  subroutine dpbd( &
       mfbe,img, &
       typdat, &
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
!_A    Realisation d'ecritures de donnees sur des points frontieres.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    use mod_dpbdn
    use mod_dpbdc
    use mod_dpbdb
    use mod_dpbdi
    implicit none
    integer          ::        img,         l,        lb,      mfbe, mnc(ip43)
    integer          ::  mnr(ip44),ncbd(ip41),ncin(ip41)
    double precision :: nxn(ip42),nyn(ip42),nzn(ip42),tm1(ip40),tm2(ip40)
    double precision :: tm3(ip40),tm4(ip40),tm5(ip40),tm6(ip40),xnr(ip44)
    double precision :: ynr(ip44),znr(ip44)
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  typdat
!$OMP MASTER
!
!
    l=ndlb(nfei(mfbe))
!
    if (typdat(1:2).eq.'st') then
       call dpbdi( &
            mfbe,img, &
            ncbd)
    elseif (typdat(1:5).eq.'basic') then
       call dpbdb( &
            mfbe,img, &
            ncbd,ncin)
    elseif (typdat(1:4).eq.'norm') then
       call dpbdn( &
            mfbe,img, &
            ncbd,nxn,nyn,nzn)
    elseif (typdat(1:4).eq.'coin') then
       lb=ndcc(nfei(mfbe))
       call dpbdc( &
            mfbe,img, &
            mnc, &
            tm1,tm2,tm3)
    elseif (typdat(1:7).eq.'noncoin') then
!            call dpbdr( &
!                 mfbe,img, &
!                 mnr,xnr,ynr,znr, &
!                 tm1,tm2,tm3, &
!                 tm4,tm5,tm6)
    else
       stop 'Type de sortie frontiere inconnu'
    endif
!
!$OMP END MASTER
    return
  end subroutine dpbd
end module mod_dpbd
