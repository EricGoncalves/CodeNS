module mod_sch_duup
  implicit none
contains
  subroutine sch_duup( &
       sn,vol,t, &
       v,ptdual,vdual,vdual1, &
       mu,mut,dist,tprod, &
       ncyc,img, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DA  DATE_C : octobre 2001 - Eric Goncalves - SINUMEF
!
!     ACT
!_A    Remise a jour des valeurs  v et de ptdual
!
!----------------------------------------------------------------
! v             champ a l'instant courant
! vdual         champ a l'instant n
! vdual1        champ a l'instant n-1
! ptdual        Partie de la derivee temporelle a l'instant courant
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use modeleturb
    use schemanum
    use sortiefichier
    use chainecarac
    use mod_met_klrmut
    use mod_met_klmut
    use mod_met_samut
    use mod_met_klsmut
    use mod_at_cutke
    use mod_met_kocmut
    use mod_met_kemutr
    use mod_met_kemutm
    use mod_met_komut
    use mod_met_kemut
    use mod_met_komutr
    use mod_met_klnmut
    implicit none
    integer          ::    i,  i1,  i2,i2m1, img
    integer          :: ind1,ind2,   j,  j1,  j2
    integer          :: j2m1,   k,  k1,  k2,k2m1
    integer          ::    l,  lm,   m, n0c,  nc
    integer          :: ncyc, nid,nijd, njd
    double precision ::       cmui1(ip21),      cmui2(ip21),      cmuj1(ip21),      cmuj2(ip21),      cmuk1(ip21)
    double precision ::       cmuk2(ip21),       dist(ip12),         mu(ip12),        mut(ip12),ptdual(ip11,ip60)
    double precision ::     sn(ip31*ndir),          t(ip00),      tprod(ip00),     v(ip11,ip60), vdual(ip11,ip60)
    double precision :: vdual1(ip11,ip60),        vol(ip11)
    double precision,allocatable :: dvxx(:),dvxy(:),dvxz(:),dvyx(:),dvyy(:)
    double precision,allocatable :: dvyz(:),dvzx(:),dvzy(:),dvzz(:)
!
!-----------------------------------------------------------------------
!
!
!


    ALLOCATE(dvxx(ip00),dvxy(ip00),dvxz(ip00),dvyx(ip00),dvyy(ip00),dvyz(ip00), &
         dvzx(ip00),dvzy(ip00),dvzz(ip00))

    do l = 1,lzx
       lm=l+(img-1)*lz
       n0c=npc(lm)
       i1=ii1(lm)
       i2=ii2(lm)
       j1=jj1(lm)
       j2=jj2(lm)
       k1=kk1(lm)
       k2=kk2(lm)
       nid=id2(lm)-id1(lm)+1
       njd=jd2(lm)-jd1(lm)+1
       nijd=nid*njd
       i2m1=i2-1
       j2m1=j2-1
       k2m1=k2-1
       ind1=indc(i1  ,j1  ,k1  )
       ind2=indc(i2m1,j2m1,k2m1)
!
!     Stockage de  vdual dans vdual1 et v dans vdual
!     Calcul de ptdual
       do m=ind1,ind2
          nc=m+n0c
          vdual1(nc,1) = vdual(nc,1)
          vdual1(nc,2) = vdual(nc,2)
          vdual1(nc,3) = vdual(nc,3)
          vdual1(nc,4) = vdual(nc,4)
          vdual1(nc,5) = vdual(nc,5)
!
          vdual(nc,1) = v(nc,1)
          vdual(nc,2) = v(nc,2)
          vdual(nc,3) = v(nc,3)
          vdual(nc,4) = v(nc,4)
          vdual(nc,5) = v(nc,5)
!
          ptdual(nc,1) = -2.*vdual(nc,1) + 0.5*vdual1(nc,1)
          ptdual(nc,2) = -2.*vdual(nc,2) + 0.5*vdual1(nc,2)
          ptdual(nc,3) = -2.*vdual(nc,3) + 0.5*vdual1(nc,3)
          ptdual(nc,4) = -2.*vdual(nc,4) + 0.5*vdual1(nc,4)
          ptdual(nc,5) = -2.*vdual(nc,5) + 0.5*vdual1(nc,5)
       enddo
!
       if(kdualns.gt.0) then
          call at_cutke(lm,v)
!
          do m=ind1,ind2
             nc=m+n0c
             vdual1(nc,6) = vdual(nc,6)
             vdual1(nc,7) = vdual(nc,7)
             vdual(nc,6) = v(nc,6)
             vdual(nc,7) = v(nc,7)
             ptdual(nc,6) = -2.*vdual(nc,6) + 0.5*vdual1(nc,6)
             ptdual(nc,7) = -2.*vdual(nc,7) + 0.5*vdual1(nc,7)
          enddo
!
          if(kcmut.eq.1) then
!        Jones Launder ou Launder Sharma
!        modele de base et avec correction SST
             if((equatt(4:4).eq.' ').or.(equatt(4:4).eq.'S').or. &
                  (equatt(4:4).eq.'C').or.(equatt(4:4).eq.'L')) then
!            modele de base ou avec correction SST
                call met_kemut( &
                     lm, &
                     sn,vol,t, &
                     dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                     dist,v,mu,mut, &
                     cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
             else if(equatt(4:4).eq.'R') then
!            modele de Jones Launder realisable
                call met_kemutr( &
                     lm,ncyc, &
                     sn,vol,t, &
                     dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                     v,mu,mut, &
                     cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
             else if(equatt(4:4).eq.'M') then
!            modele de Jones Launder modifie avec C_mu variable
                call met_kemutm( &
                     lm, &
                     sn,vol,t, &
                     dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                     v,mu,mut, &
                     cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
             endif
!
          else if(kcmut.eq.2) then
!
!           Modele de Spalart-Allmaras
!
             call met_samut(lm,v,mu,mut)
!
          else if(kcmut.eq.3) then
!           Modele k-l de Smith
!
             if((equatt(4:4).eq.' ').or.(equatt(4:4).eq.'L')) then
!            modele k-l de base et SA
                call met_klmut( &
                     lm, &
                     v,mu,mut,dist)
             else if(equatt(4:4).eq.'S') then
!            modele k-l avec correction SST
                call met_klsmut( &
                     lm, &
                     sn,vol,t, &
                     dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                     dist,v,mu,mut, &
                     cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
             else if(equatt(4:4).eq.'R') then
!            modele k-l avec correction de réalisabilité
                call met_klrmut( &
                     lm, &
                     sn,vol,t, &
                     dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                     dist,v,mu,mut, &
                     cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
             else if(equatt(4:4).eq.'N') then
!            modele k-l avec correction de non equilibre
                call met_klnmut( &
                     lm, &
                     sn, &
                     vol,v,mu,mut,dist, &
                     dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                     tprod,t, &
                     cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
             else
                write(imp,'(/,''!!!sch_duup: modele k-l de Smith '',''non prevu'')')
                stop
             endif
!
          else if(kcmut.eq.4) then
!com        Modele k-omega de Menter avec realisabilite de Durbin
             if((equatt(4:4).eq.'R').or. &
                  (equatt(5:5).eq.'R')) then
                call met_komutr( &
                     lm, &
                     sn,vol,t, &
                     dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                     dist,v,mu,mut, &
                     cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
             else
!            Modeles k-omega de Wilcox et Menter
                call met_komut( &
                     lm, &
                     sn,vol,t, &
                     dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                     dist,v,mu,mut, &
                     cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
             endif
!
          else if(kcmut.eq.5) then
!           modele de Chien
!            call met_chmut(lm,v,mu,mut,dist,mnpar,utau)
!
             write(imp,'(/,''!!!sch_duup: modele k-e de Chien '',''non prevu'')')
             stop
!         else if(kcmut.eq.6) then
!
!           Jones Launder ou Launder Sharma avec correction pour utiliser
!           la fonction f_mu de Smith dans les regions externes et les sillages
!
!            call met_mutke2(
!     &           l,ncyc,
!     &           v,mu,mut,dist,mnpar,ncin,
!     &           txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z,
!     &           qcz000)
!
          else if(kcmut.eq.7) then
!        RNG avec k-l a la paroi
!            call met_rngmut(lm,v,mu,mut,dist,fracmod)
             write(imp,'(/,''!!!sch_duup: modele k-e RNG '',''non prevu'')')
             stop
!
          else if(kcmut.eq.8) then
!c         modele k-omega de Wilcox compressible
             if((equatt(4:4).eq.' ').or.(equatt(4:4).eq.'S')) then
                call met_kocmut( &
                     lm, &
                     sn,vol,t, &
                     dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                     dist,v,mu,mut, &
                     cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!            call met_kokmut(
!     &           lm,
!     &           sn,vol,t,
!     &           dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz,
!     &           dist,v,mu,mut,
!     &           cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
             else if(equatt(4:4).eq.'R') then
!            call met_kokmutr( &
!                 lm, &
!                 sn,vol,t, &
!                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
!                 v,mu,mut, &
!                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
             else
                write(imp,'(/,''!!!sch_duup: kcmut non prevu'')')
                stop
             endif
          endif      !fin test kcmut
       endif       !fin test kdualns
    enddo        !fin boucle domaine

    DEALLOCATE(dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz)

    return
  contains
    function    indc(i,j,k)
      implicit none
      integer          ::    i,indc,   j,   k
      indc=1+(i-id1(lm))+(j-jd1(lm))*nid+(k-kd1(lm))*nijd
    end function indc
  end subroutine sch_duup
end module mod_sch_duup
