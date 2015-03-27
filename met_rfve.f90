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
    integer          ::          m,        mb,        mf,       mfb,        mt
    integer          :: ncbd(ip41),ncin(ip41),        nd,        ni
    double precision :: t(ip11,ip60)
!
!-----------------------------------------------------------------------
!
!
    do mf=1,nbd
       mfb=lbd(mf)
       mt=mmb(mfb)
!!$OMP SIMD
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
  end subroutine met_rfve
end module mod_met_rfve
