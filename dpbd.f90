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
implicit none
integer :: mfbe
integer :: img
integer :: ncbd
integer :: ncin
integer :: mnc
integer :: mnr
double precision :: xnr
double precision :: ynr
double precision :: znr
double precision :: tm1
double precision :: tm2
double precision :: tm3
double precision :: tm4
double precision :: tm5
double precision :: tm6
integer :: l
integer :: lb
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  typdat
      real nxn,nyn,nzn
!
      dimension ncbd(ip41),ncin(ip41)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
      dimension mnc(ip43)
      dimension xnr(ip44),ynr(ip44),znr(ip44),mnr(ip44)
      dimension tm1(ip40),tm2(ip40),tm3(ip40),tm4(ip40), &
                tm5(ip40),tm6(ip40)
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
      return
      end
