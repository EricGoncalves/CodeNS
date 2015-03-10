module mod_met_rbse
implicit none
contains
      subroutine met_rbse( &
                 s1x,s1y,s1z,s2x,s2y,s2z, &
                 ncbd,ncin)
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
double precision :: s1x
double precision :: s1y
double precision :: s1z
double precision :: s2x
double precision :: s2y
double precision :: s2z
integer :: ncbd
integer :: ncin
integer :: m
integer :: mb
integer :: mf
integer :: mfb
integer :: mt
integer :: nd
integer :: ndm
!
!-----------------------------------------------------------------------
!
      dimension s1x(ip12),s1y (ip12),s1z(ip12), &
                s2x(ip12),s2y (ip12),s2z(ip12)
      dimension ncin(ip41),ncbd(ip41)
!
      do mf=1,nbd
!
         mfb=lbd(mf)
         mt=mmb(mfb)
!
!DEC$ IVDEP
         do m=1,mt
            mb=mpb(mfb)+m
            nd=ncbd(mb)
            ndm=ncin(mb)
!
!     definition des variables aux bords (centre des facettes frontieres)
!
            s1x(nd) =s1x(ndm)
            s1y(nd) =s1y(ndm)
            s1z(nd) =s1z(ndm)
            s2x(nd) =s2x(ndm)
            s2y(nd) =s2y(ndm)
            s2z(nd) =s2z(ndm)
!
         enddo
!
      enddo
!
      return
      end
end module
