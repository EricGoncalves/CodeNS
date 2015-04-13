module mod_met_rbve
  implicit none
contains
  subroutine met_rbve(t,ncin,ncbd)
!
!***********************************************************************
!
!     ACT  extrapolation des grandeurs turbulentes aux bords
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    implicit none
    integer          ::          m,        mf,       mfb,        ml,        mt
    integer          ::          n,ncbd(ip41),ncin(ip41),        ni
    double precision :: t(ip11,ip60)
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!
!     definition des variables aux bords (centre des facettes frontieres)
!
    do mf=1,nbd
       mfb=lbd(mf)
       mt=mmb(mfb)
       do m=1,mt
          ml=mpb(mfb)+m
          n=ncbd(ml)
          ni=ncin(ml)
          t(n,6) = t(ni,6)
          t(n,7) = t(ni,7)
       enddo
    enddo
!
!$OMP END MASTER
    return
  end subroutine met_rbve
end module mod_met_rbve
