module mod_sch_duin
  implicit none
contains
  subroutine sch_duin( &
       v,img,ptdual,vdual,vdual1,vdual2)
!
!***********************************************************************
!
!_DA           octobre 2001 - Eric Goncalves - SINUMEF
!
!     ACT
!_A    Initialisation de ptdual sur grille fine.
!
!----------------------------------------------------------------
! v             champ a l'instant courant
! vdual         champ a l'instant n
! vdual1        champ a l'instant n-1
! vdual2        champ a l'instant n-2
! ptdual        Partie de la derivee temporelle a l'instant courant
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use schemanum
    implicit none
    integer          ::    i,  i1,  i2,i2m1, img
    integer          :: ind1,ind2,   j,  j1,  j2
    integer          :: j2m1,   k,  k1,  k2,k2m1
    integer          ::    l,  lm,   m, n0c,  nc
    integer          ::  nid,nijd, njd
    double precision ::                c1,ptdual(ip11,ip60),     v(ip11,ip60), vdual(ip11,ip60),vdual1(ip11,ip60)
    double precision :: vdual2(ip11,ip60)
!
!-----------------------------------------------------------------------
!
!

    c1=1./3.
!
    do l = 1,lzx
       lm = l + (img-1)*lz
       n0c=npc(lm)
       i1=ii1(lm)
       i2=ii2(lm)
       j1=jj1(lm)
       j2=jj2(lm)
       k1=kk1(lm)
       k2=kk2(lm)
       nid = id2(lm)-id1(lm)+1
       njd = jd2(lm)-jd1(lm)+1
       nijd = nid*njd
       i2m1=i2-1
       j2m1=j2-1
       k2m1=k2-1
       ind1 = indc(i1  ,j1  ,k1  )
       ind2 = indc(i2m1,j2m1,k2m1)
!
       do m=ind1,ind2
          nc=m+n0c
!
! Initialisation a l'ORDRE 1 de la derivee en temps physique
!
          ptdual(nc,1) = -vdual(nc,1)
          ptdual(nc,2) = -vdual(nc,2)
          ptdual(nc,3) = -vdual(nc,3)
          ptdual(nc,4) = -vdual(nc,4)
          ptdual(nc,5) = -vdual(nc,5)
          ptdual(nc,6) = -vdual(nc,6)
          ptdual(nc,7) = -vdual(nc,7)
!
! Initialisation a l'ORDRE 2 de la derivee en temps physique
!
!        ptdual(nc,1) = -2.*vdual(nc,1) + 0.5*vdual1(nc,1)
!        ptdual(nc,2) = -2.*vdual(nc,2) + 0.5*vdual1(nc,2)
!        ptdual(nc,3) = -2.*vdual(nc,3) + 0.5*vdual1(nc,3)
!        ptdual(nc,4) = -2.*vdual(nc,4) + 0.5*vdual1(nc,4)
!        ptdual(nc,5) = -2.*vdual(nc,5) + 0.5*vdual1(nc,5)
!        ptdual(nc,6) = -2.*vdual(nc,6) + 0.5*vdual1(nc,6)
!        ptdual(nc,7) = -2.*vdual(nc,7) + 0.5*vdual1(nc,7)
!
! Initialisation a l'ORDRE 3 de la derivee en temps physique
!
!         v(nc,1) = vdual(nc,1) &
!        + (11.0*vdual(nc,1) - 18.0*vdual1(nc,1) + 9.0*vdual2(nc,1)- 2.0*vdual3(nc,1))/6.
!         v(nc,2) = vdual(nc,2) &
!        + (11.0*vdual(nc,2) - 18.0*vdual1(nc,2) + 9.0*vdual2(nc,2)- 2.0*vdual3(nc,2))/6.
!         v(nc,3) = vdual(nc,3) &
!        + (11.0*vdual(nc,3) - 18.0*vdual1(nc,3) + 9.0*vdual2(nc,3)- 2.0*vdual3(nc,3))/6.
!         v(nc,4) = vdual(nc,4) &
!        + (11.0*vdual(nc,4) - 18.0*vdual1(nc,4) + 9.0*vdual2(nc,4)- 2.0*vdual3(nc,4))/6.
!         v(nc,5) = vdual(nc,5) &
!        + (11.0*vdual(nc,5) - 18.0*vdual1(nc,5) + 9.0*vdual2(nc,5)- 2.0*vdual3(nc,5))/6.
!         v(nc,6) = vdual(nc,6) &
!        + (11.0*vdual(nc,6) - 18.0*vdual1(nc,6) + 9.0*vdual2(nc,6)- 2.0*vdual3(nc,6))/6.
!         v(nc,5) = vdual(nc,5) &
!        + (11.0*vdual(nc,7) - 18.0*vdual1(nc,7) + 9.0*vdual2(nc,7)- 2.0*vdual3(nc,7))/6.
!
!         ptdual(nc,1) = -3.*vdual(nc,1) + 1.5*vdual1(nc,1) - c1*vdual2(nc,1)
!         ptdual(nc,2) = -3.*vdual(nc,2) + 1.5*vdual1(nc,2) - c1*vdual2(nc,2)
!         ptdual(nc,3) = -3.*vdual(nc,3) + 1.5*vdual1(nc,3) - c1*vdual2(nc,3)
!         ptdual(nc,4) = -3.*vdual(nc,4) + 1.5*vdual1(nc,4) - c1*vdual2(nc,4)
!         ptdual(nc,5) = -3.*vdual(nc,5) + 1.5*vdual1(nc,5) - c1*vdual2(nc,5)
!         ptdual(nc,6) = -3.*vdual(nc,6) + 1.5*vdual1(nc,6) - c1*vdual2(nc,6)
!         ptdual(nc,7) = -3.*vdual(nc,7) + 1.5*vdual1(nc,7) - c1*vdual2(nc,7)
       enddo
    enddo
!
    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
    end function indc
  end subroutine sch_duin
end module mod_sch_duin
