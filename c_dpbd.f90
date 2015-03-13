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
integer :: imot
integer :: nmot
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
integer :: img
integer :: lgr
integer :: lgrd
integer :: lmfb
integer :: lmfbd
integer :: mf
integer :: mfbe
integer :: mg
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot(nmx)
      character(len=32) ::  typdat
      double precision nxn,nyn,nzn
!
      dimension imot(nmx)
      dimension lmfb(mtb)
      dimension lgr(nobj)
!
      dimension ncbd(ip41),ncin(ip41)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
      dimension mnc(ip43)
      dimension xnr(ip44),ynr(ip44),znr(ip44),mnr(ip44)
      dimension tm1(ip40),tm2(ip40),tm3(ip40),tm4(ip40), &
                tm5(ip40),tm6(ip40)
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
      end subroutine
end module
