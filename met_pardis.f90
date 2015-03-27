module mod_met_pardis
  implicit none
contains
  subroutine met_pardis(ncin,d)
!
!***********************************************************************
!
!_DA  DATE_C : avril 1999-- AUTEUR : Eric Goncalves
!
!     ACT
!_A   Mise a zero de la dissipation a la paroi
!_A   Modele k-omega de Wilcox et Menter
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use mod_atindnor
    implicit none
    integer          ::        dm1,       dm2,       dm3,    isens3,         m
    integer          ::        m10,     m1max,     m1min,       m20,     m2max
    integer          ::      m2min,        m3,       m30,     m3max,     m3min
    integer          ::         mb,       mfb,        mt,         n,ncin(ip41)
    integer          ::         ni,        no,      nper
    double precision :: d(ip11,ip60)
!
!-----------------------------------------------------------------------
!
!
!c    boucle sur toutes les frontieres
    do no=1,mtbx
       mfb=nba(no)
       if (cl(mfb)(1:3).eq.'par') then
          mt=mmb(mfb)
!c      suppression de la dissipation artificielle pour omega pres
!c      de la paroi (sur nper couches de cellules perpendiculairement
!c      a la paroi)
!c      calculs des indices direction parallele et normale a la paroi
          call atindnor( &
               mfb, &
               m10,m20,m30, &
               m1min,m1max,m2min,m2max,m3min,m3max, &
               dm1,dm2,dm3,isens3)
!
          nper=3
          m3max=m3min+isens3*(nper-1)
          do m=1,mt
             mb=mpb(mfb)+m
             ni=ncin(mb)
             do m3=m3min,m3max,isens3
!c         boucle sur les nper points perpendiculaire a la paroi
                n=ni+(m3-m3min)*dm3
                d(n,7)=0.
             enddo
          enddo
       endif
    enddo
!
    return
  end subroutine met_pardis

end module mod_met_pardis
