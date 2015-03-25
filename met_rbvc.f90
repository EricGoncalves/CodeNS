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
  integer          ::          m,        mb,        mc,        mf,       mfb
  integer          ::  mnc(ip43),        mt,        nc,ncbd(ip41),ncin(ip41)
  integer          ::         nd,       ndm
  double precision :: t(ip11,ip60)
!
!-----------------------------------------------------------------------
!
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
  end subroutine met_rbvc
end module mod_met_rbvc
