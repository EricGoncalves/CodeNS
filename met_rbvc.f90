module mod_met_rbvc
implicit none
contains
      subroutine met_rbvc( &
                 t, &
                 ncbd,ncin,mnc)
!
!***********************************************************************
!
!     ACT
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
integer :: mnc
integer :: m
integer :: mb
integer :: mc
integer :: mf
integer :: mfb
integer :: mt
integer :: nc
integer :: nd
integer :: ndm
!
!-----------------------------------------------------------------------
!
      dimension t(ip11,ip60)
      dimension ncbd(ip41)
      dimension ncin(ip41)
      dimension mnc(ip43)
!
      do mf=1,nbd
!
         mfb=lbd(mf)
         mt=mmb(mfb)
!
!!$OMP SIMD
         do m=1,mt
            mc=mpc(mfb)+m
            nc=mnc(mc)
            mb=mpb(mfb)+m
            nd=ncbd(mb)
            ndm=ncin(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
            t(nd,6) = 0.5*( t(ndm,6)+t(nc,6) )
            t(nd,7) = 0.5*( t(ndm,7)+t(nc,7) )
!
         enddo
!
      enddo
!
      return
      end subroutine
end module
