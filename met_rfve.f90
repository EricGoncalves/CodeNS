module mod_met_rfve
implicit none
contains
      subroutine met_rfve(t,ncbd,ncin)
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
   use boundary
implicit none
double precision :: t
integer :: ncbd
integer :: ncin
integer :: m
integer :: mb
integer :: mf
integer :: mfb
integer :: mt
integer :: nd
integer :: ni
!
!-----------------------------------------------------------------------
!
      dimension t(ip11,ip60)
      dimension ncin(ip41),ncbd(ip41)
!
      do mf=1,nbd
       mfb=lbd(mf)
       mt=mmb(mfb)
!DEC$ IVDEP
       do m=1,mt
        mb=mpb(mfb)+m
        nd=ncbd(mb)
        ni=ncin(mb)
        t(nd,6)=t(ni,6)
        t(nd,7)=t(ni,7)
       enddo
      enddo
!
      return
      end subroutine
end module
