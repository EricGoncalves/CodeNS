module mod_met_num
  use mod_zvismo
  use mod_met_uttau
  implicit none
contains
  subroutine met_num( &
       ncbd,ncin,mnc,ncyc, &
       mnr,xnr,ynr,znr, &
       bceqt, &
       u,v,dt,d,mut,mu,cfke, &
       mnpar,dist,fgam,utau,topz, &
       nxn,nyn,nzn, &
       txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
       qcxts5,qcyts6,qcz000, &
       sn,vol, &
       icycle,ptdual, &
       x,y,z, &
       tn1,tn2,tn3,t,tprod,dtdx,dtdy,dtdz,tn9,tn10, &
       tp,cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
       pression,ztemp,cson)
!
!***********************************************************************
!
!_DA :  AUTEUR: Eric Goncalves / SINUMEF
!
!     ACT
!_A   Integration de 2 equations de transport de n a n+1 par schema
!_A   centre de Jameson ou decentre de Roe.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use modeleturb
    use schemanum
    use proprieteflu
    use chainecarac
    use sortiefichier
    use mod_lp2ke
    use mod_met_samut
    use mod_lp2kl
    use mod_met_chmut
    use mod_met_fludcsa
    use mod_met_roe2oh
    use mod_met_smsa
    use mod_lp2kw
    use mod_met_smklsas
    use mod_sch_turb2
    use mod_met_smdes
    use mod_impli2_eqt_3d
    use mod_sch_turb
    use mod_met_komutr
    use mod_met_roe2o
    use mod_met_kocmut
    use mod_met_bord
    use mod_lpker
    use mod_at_cutke
    use mod_met_roe
    use mod_met_mutke2
    use mod_met_fludmt
    use mod_lpke
    use mod_met_klrmut
    use mod_met_klsmut
    use mod_atctranske
    use mod_lpkl
    use mod_met_dual2
    use mod_lp2sa
    use mod_lpkomega
    use mod_met_gradtr
    use mod_met_fludc
    use mod_met_smch
    use mod_met_smsasas
    use mod_met_kemutr
    use mod_met_fludko
    use mod_met_smmtr
    use mod_impli2_eqt
    use mod_met_kemutm
    use mod_met_smkor
    use mod_met_rbsc
    use mod_met_rbse
    use mod_met_parko
    use mod_lpsa
    use mod_met_smmt
    use mod_met_klmut
    use mod_met_komut
    use mod_met_prod
    use mod_met_dual
    use mod_met_smke
    use mod_met_rbsr
    use mod_lpkomegar
    use mod_met_kemut
    use mod_met_smkl
    use mod_met_smko
    implicit none
    integer          ::      icycle,          l,     lgsnlt,         mf,        mfc
    integer          ::         mfr,  mnc(ip43),mnpar(ip12),  mnr(ip44), ncbd(ip41)
    integer          ::  ncin(ip41),       ncyc,       npsn
    double precision ::  bceqt(ip41,neqt),       cfke(ip13),      cmui1(ip21),      cmui2(ip21),      cmuj1(ip21)
    double precision ::       cmuj2(ip21),      cmuk1(ip21),      cmuk2(ip21),       cson(ip11),     d(ip11,ip60)
    double precision ::        dist(ip12),         dt(ip11),       dtdx(ip00),       dtdy(ip00),       dtdz(ip00)
    double precision ::        fgam(ip42),         mu(ip12),        mut(ip12),        nxn(ip42),        nyn(ip42)
    double precision ::         nzn(ip42),   pression(ip11),ptdual(ip11,ip60),     qcxts5(ip12),     qcyts6(ip12)
    double precision ::      qcz000(ip12),    sn(ip31*ndir),          t(ip00),        tn1(ip00),       tn10(ip00)
    double precision ::         tn2(ip00),        tn3(ip00),        tn9(ip00),       topz(ip11),         tp(ip40)
    double precision ::       tprod(ip00),     txxf5x(ip12),     txyf5y(ip12),     txzf5z(ip12),     tyyf6x(ip12)
    double precision ::      tyzf6y(ip12),     tzzf6z(ip12),     u(ip11,ip60),       utau(ip42),     v(ip11,ip60)
    double precision ::         vol(ip11),          x(ip21),        xnr(ip44),          y(ip21),        ynr(ip44)
    double precision ::           z(ip21),        znr(ip44),      ztemp(ip11)
    double precision,allocatable :: coefe(:,:),   dvxx(:),   dvxy(:),   dvxz(:),   dvyx(:)
    double precision,allocatable ::    dvyy(:),   dvyz(:),   dvzx(:),   dvzy(:),   dvzz(:)
    double precision,allocatable :: fracmod(:)
!$OMP MASTER
    allocate(coefe(ndir,ip00))
!
!-----------------------------------------------------------------------
!
!
    ALLOCATE(dvxx(ip00),dvxy(ip00),dvxz(ip00),dvyx(ip00),dvyy(ip00),dvyz(ip00), &
         dvzx(ip00),dvzy(ip00),dvzz(ip00),fracmod(ip12))
!
!     choix des appels : a initialiser en fonction du modele
!     Fait dans "atlecdon.f"
!     kcutke =1
!     kfludis=1
!     ksecmb =1
!     kcmut  =1
!
    if((kutau.eq.1).and.(ncyc.ne.1)) then
!
!          modele de Chien - modele k-omega bas Reynolds de Wilcox
!          et modeles k-omega de Wilcox ou Menter avec rugosite
!          il faut "utau" a partir du tenseur des contraintes
!          visqueuses et turbulentes.
!
       call met_uttau( &
            ncbd, &
            nxn,nyn,nzn, &
            txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
            v,utau)
!
    endif
!
!-----------application des condition aux limites----------------------------------------
!
    call met_bord( &
         0, &
         bceqt, &
         mnc,ncin,mnr,xnr,ynr,znr,ncbd, &
         v,utau,mu)
!
    if((kparoi.eq.1).and.(lparoi.lt.1)) then
!
!              condition limite sur omega en 6*nu/y**2 pour Wilcox et Menter
!              traitement de toutes les surfaces
!
       call met_parko( &
            0,ncin,ncbd, &
            v,mut,mu,dist,mnpar)
    endif
!
! --------------------------------------------------------------------------
    do l=1,lzx
!
!          calcul de la viscosite moleculaire
       call zvismo(l,mu,v,ztemp)

!           ----------------------------------------------------------
!com        teq_gradv --> grad(v) aux points interieurs au domaine
!
!           | dvxx,dvxy,dvxz |    | du/dx du/dy du/dz |
!           | dvyx,dvyy,dvyz | => | dv/dx dv/dy dv/dz |
!           | dvzx,dvzy,dvzz |    | dw/dx dw/dy dw/dz |
!
       call teq_gradv( &
            l, &
            sn, &
            vol,v, &
            t, &
            dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
            cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!           ----------------------------------------------------------
!com        met_prod --> calcul de v x production de k
!
       call met_prod( &
            l,ncyc, &
            mut,v, &
            dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
            tprod)
!           ----------------------------------------------------------
!com        application de la condition aux limites lois de paroi
!
       if(lparoi.eq.1) then
          if((equatt(1:3).eq.'2JL').or.(equatt(1:3).eq.'2LS')) then
!              Modele k-epsilon de Jones Launder
!              production de k et epsilon dans cellule adjacente aux parois
             if((equatt(4:4).eq.' ').or.(equatt(4:4).eq.'S').or. &
                  (equatt(4:4).eq.'C').or.(equatt(4:4).eq.'L')) then
                call lpke( &
                     v,mu,mut,dist, &
                     nxn,nyn,nzn, &
                     ncin,ncbd,l, &
                     mnpar,fgam,ncyc, &
                     tprod,tp, &
                     ztemp)
!
             elseif (equatt(4:4).eq.'R') then
                call lpker( &
                     v,mu,mut,dist, &
                     nxn,nyn,nzn, &
                     ncin,ncbd,l, &
                     mnpar,fgam,tprod, &
                     dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                     ztemp)
             endif
!
          else if(equatt(1:3).eq.'2Sm') then
!              Modele k-l de Smith
!              production de k et l imposee dans cellule adjacente aux parois
             call lpkl( &
                  v,mu,mut,dist, &
                  nxn,nyn,nzn, &
                  ncin,ncbd,l, &
                  mnpar,fgam, &
                  tprod,ncyc,tp, &
                  ztemp)
!
          else if((equatt(1:3).eq.'2WL').or.(equatt(1:3).eq.'2KO') &
               .or.(equatt(1:3).eq.'2MT')) then
!             Modele k-omega de Wilcox et Menter
!             production de k et omega dans cellule adjacente aux parois
             if((equatt(4:4).eq.' ').or.(equatt(4:4).eq.'S')) then
                call lpkomega( &
                     v,mu,mut,dist, &
                     nxn,nyn,nzn, &
                     ncin,ncbd,l, &
                     mnpar,fgam, &
                     tprod,ncyc,tp, &
                     ztemp)
!
             elseif(equatt(4:4).eq.'G') then  !earsm gatski
!                 call lpkomegag( &
!                 v,mu,mut,dist, &
!                 nxn,nyn,nzn, &
!                 ncin,ncbd,l, &
!                 mnpar,fgam,tprod,tp, &
!                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz)
!
             elseif((equatt(4:4).eq.'R').or. &
                  (equatt(5:5).eq.'R')) then    !realisabilite durbin
                call lpkomegar( &
                     v,mu,mut,dist, &
                     nxn,nyn,nzn, &
                     ncin,ncbd,l, &
                     mnpar,fgam,tprod,tp, &
                     dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                     ztemp)
             endif
!
          else if(equatt(1:3).eq.'1SA') then
!              Modele de Spalart-Allmaras
!              nutilde imposee dans cellule adjacente aux parois
             call lpsa( &
                  v,mu,mut,dist, &
                  nxn,nyn,nzn, &
                  ncin,ncbd,l, &
                  mnpar,fgam,ncyc,tp, &
                  ztemp)
          endif
!
       elseif(lparoi.eq.2) then
          if(ncyc.ne.1) then
             call met_uttau( &
                  ncbd, &
                  nxn,nyn,nzn, &
                  txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                  v,utau)
          endif
          if(equatt(1:3).eq.'2Sm') then
!              Modele k-l de Smith
             call lp2kl( &
                  v,mu,mut,dist, &
                  nxn,nyn,nzn, &
                  ncin,ncbd,l, &
                  vol,sn,ncyc, &
                  mnpar,fgam,  &
                  tprod,tp,utau,topz, &
                  tn1,tn2,tn3, &
                  pression,ztemp)
          elseif((equatt(1:3).eq.'2JL').or.(equatt(1:3).eq.'2LS')) then
!              Modele k-eps de Jones-Launder ou de Launder-Sharma
             call lp2ke( &
                  v,mu,mut,dist, &
                  nxn,nyn,nzn, &
                  ncin,ncbd,l, &
                  vol,sn,ncyc, &
                  mnpar,fgam,  &
                  tprod,tp,utau, &
                  tn1,tn2,tn3, &
                  pression,ztemp)
          elseif((equatt(1:3).eq.'2WL').or.(equatt(1:3).eq.'2KO') &
               .or.(equatt(1:3).eq.'2MT')) then
!              Modele k-omega de Wilcox ou Menter
             call lp2kw( &
                  v,mu,mut,dist, &
                  nxn,nyn,nzn, &
                  ncin,ncbd,l, &
                  vol,sn,ncyc, &
                  mnpar,fgam,  &
                  tprod,tp,utau, &
                  tn1,tn2,tn3, &
                  pression,ztemp)
          elseif(equatt(1:3).eq.'1SA') then
!              Modele de Spalart-Allmaras
             call lp2sa( &
                  v,mu,mut,dist, &
                  nxn,nyn,nzn, &
                  ncin,ncbd,l, &
                  vol,sn,ncyc, &
                  mnpar,fgam,  &
                  tp,utau,     &
                  tn1,tn2,tn3, &
                  pression,ztemp)
          endif
       endif
!           ----------------------------------------------------------
!com        grad(k)=>txxf5x,txyf5y,txzf5z points interieurs au domaine
!com        grad(e)=>tyyf6x,tyzf6y,tzzf6z
!
       call met_gradtr( &
            l, &
            sn, &
            vol,v,mu,mut, &
            t,dtdx,dtdy,dtdz, &
            txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
            cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!           -----------------------------------------------------------
!
!com        calcul terme source des equations de transport
!
       if(ksecmb.le.1) then
!
!              modele de Jones Launder
!
          call met_smke( &
               l, &
               sn, &
               vol,v,mu,mut,tprod,cfke, &
               dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
               t,dtdx,dtdy,dtdz,tn9,tn10, &
               qcxts5,qcyts6, &
               tn1,tn2,tn3, &
               cson, &
               cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
       else if(ksecmb.eq.2) then
!
          if(equatt(1:4).eq.'1SA ') then
!
!              modele de Spalart-Allmaras
!
             call met_smsa( &
                  l, &
                  sn, &
                  vol,v,mu,dist, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  txxf5x,txyf5y,txzf5z,cfke, &
                  t,dtdx,dtdy,dtdz,tn9, &
                  qcxts5,qcyts6, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
          elseif(equatt(1:4).eq.'1SAL') then
!
!              modele de Spalart-Allmaras SAS
!
             call met_smsasas( &
                  l, &
                  sn, &
                  vol,v,mu,dist, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  txxf5x,txyf5y,txzf5z,cfke,  &
                  t,dtdx,dtdy,dtdz,tn9, &
                  qcxts5,qcyts6, &
                  tn1,tn2,tn3,tn10, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
          elseif(equatt(1:4).eq.'1SAD') then
!
!              DES de Spalart
!
             call met_smdes( &
                  l,x,y,z,tn9, &
                  sn, &
                  vol,v,mu,dist, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  txxf5x,txyf5y,txzf5z,cfke, &
                  t,dtdx,dtdy,dtdz,tn10, &
                  qcxts5,qcyts6, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          endif
!
       else if(ksecmb.eq.3) then
!
!              modele k-l de Smith
!              modele k-l de Smith Scale-Adaptative
!
          if(equatt(4:4).eq.'L') then
             call met_smklsas( &
                  l, &
                  sn, &
                  vol,v,mu,mut,dist,mnpar, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                  tprod,cfke, &
                  t,dtdx,dtdy,dtdz,tn10, &
                  qcxts5,qcyts6, &
                  tn1,tn2,tn3,tn9, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          else
             call met_smkl( &
                  l, &
                  sn, &
                  vol,v,mu,mut,dist,mnpar, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                  tprod,cfke, &
                  t,dtdx,dtdy,dtdz,tn9, &
                  qcxts5,qcyts6, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          endif
!
       else if(ksecmb.eq.4) then
!
!              modeles k-omega de Wilcox, Menter, Menter SST
!              modeles k-omega de Menter + realisabilite de Durbin
!
          if((equatt(4:4).eq.' ').or.(equatt(4:4).eq.'S')) then
             call met_smmt( &
                  l,ncyc, &
                  v,mu,mut,dist,mnpar,ncin, &
                  txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                  tprod,cfke, &
                  tn9,qcz000, &
                  qcxts5,qcyts6)
!
          elseif(equatt(4:4).eq.'R') then
             call met_smmtr( &
                  l, &
                  v,mu,mut,dist,mnpar,ncin, &
                  txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  tprod,cfke,fracmod, &
                  qcxts5,qcyts6)
          endif
!
       else if(ksecmb.eq.5) then
!
!              modele de Chien
!
          call met_smch( &
               l, &
               v,mu,dist,mnpar,utau, &
               tprod, &
               qcxts5,qcyts6)
!
       else if(ksecmb.eq.6) then
!
!             modele RNG de Yakhot, Orzag, Thangam, Gatski et Speziale
!
!              call met_smrng( &
!                 l,ncyc, &
!                 sn,vol, &
!                 v,mu,mut,dist,mnpar,ncin, &
!                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
!                 txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
!                 tprod,cfke, &
!                 t,dtdx,dtdy,dtdz,tn9, &
!                 fracmod,qcxts5,qcyts6, &
!                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!
       else if(ksecmb.eq.8) then
!
!             modele k-omega de Wilcox compressible
!
          if((equatt(4:4).eq.' ').or.(equatt(4:4).eq.'S')) then
             call met_smko( &
                  l,v, &
                  txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                  tprod,cfke, &
                  qcxts5,qcyts6,cson)
!
          elseif(equatt(4:4).eq.'R') then  !Kok realisable
             call met_smkor( &
                  l,v, &
                  txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  tprod,cfke, &
                  qcxts5,qcyts6)
          endif
!
       else
          write(imp,'(/,''!!!met_num: modele non prevu'')')
          stop
       endif
!
! -----------------------------------------------------------------
!        calcul des densites de flux dissipatifs
!-------------------------------------------------------------------
!           (txxf5x,txyf5y,txzf5z)<=grad(k) devient (mu+mut/.)grad(k)
!           (tyyf6x,tyzf6y,tzzf6z)<=grad(e) devient (mu+mut/.)grad(e)
!
       if(kfludis.eq.1) then !sigma_k et sigma_e constants
          call met_fludc( &
               l, &
               v,mu,mut, &
               t,dtdx,dtdy,dtdz, &
               txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z)
!
       else if(kfludis.eq.2) then !modele Spalart Allmaras
          call met_fludcsa( &
               l, &
               v,mu, &
               txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z)
!
       else if(kfludis.eq.4) then ! modele de Menter avec fonction de raccordement
          call met_fludmt( &
               l, &
               v,mu,mut, &
               qcz000, &
               txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z)
!
       else if(kfludis.eq.5) then !modele kw de Kok et Wilcox compressible
          call met_fludko( &
               l,mu,mut, &
               txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z)
!
       else if(kfludis.eq.6) then ! modele RNG de Yakhot et al
!              call met_fludrng( &
!                 l, &
!                 v,mu,mut, &
!                 fracmod, &
!                 txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z)

       else
          write(imp,'(/,''!!!met_num: kfludis non prevu'')')
          stop
       endif
!
!           fin de boucle sur les domaines
    enddo
!
!-------------------------------------------------------------------------------
!        continuite des densites de flux dissipatifs entre les domaines a tn
!-------------------------------------------------------------------------------
!
!        frontieres a normales stockees
!$OMP SIMD
    do mf=1,mtnx
       lbd(mf)=nfbn(mf)
    enddo
    nbd=mtnx
    call met_rbse( &
         txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
         ncbd,ncin)

!       frontieres autres
!$OMP SIMD
    do mf=1,mtax
       lbd(mf)=nfba(mf)
    enddo
    nbd=mtax
    call met_rbse( &
         txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
         ncbd,ncin)

!        frontieres coincidentes
!$OMP SIMD
    do mfc=1,mtcx
       lbd(mfc)=nfbc(mfc)
    enddo
    nbd=mtcx
    call met_rbsc( &
         txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
         ncbd,ncin,mnc)

!        frontiere recouverte
    do mfr=1,mtrx
       lbd(mfr)=nfbr(mfr)
    enddo
    nbd=mtrx
    call met_rbsr( &
         ncbd,mnr,xnr,ynr,znr, &
         txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
         ncin)
!
!********************************************************************
!         application du schema numerique
!********************************************************************
!
    do l=1,lzx
       npsn  =ndir*npfb(l)+1
       lgsnlt=nnn(l)
!
!-------calcul de la dissipation du schema de  Roe------------------
!
       if(kditur.eq.2) then
          if(klroe.eq.1) then
!            schema ordre 1
             call met_roe( &
                  l,v,d, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  vol)
          elseif(klroe.eq.2) then
!            schema ordre 2 - sans correction de Harten
             if(abs(epsroe).le.tiny(1.)) then
                call met_roe2o( &
                     l,v,d, &
                     sn(npsn),lgsnlt, &
                     vol,dvxx,dvxy)
             else
!            schema ordre 2 - avec correction de Harten
!            vitesse du son donne par loi d'etat
                call met_roe2oh( &
                     l,v,d, &
                     sn(npsn),lgsnlt, &
                     vol,dvxx,dvxy)
             endif
          endif
       endif
!
!------------------------------------------------------------------------
!           increment explicite
!-----------------------------------------------------------------------
!
       if((kditur.eq.1).or.(kditur.eq.2)) then
          if((ischema.eq.2).or.(ischema.eq.4).or.(ischema.eq.6).or. &
               (ischema.eq.8).or.(ischema.eq.11).or.(ischema.eq.13)) then
!             schemas ponderes de Jameson et de Roe
             call sch_turb2( &
                  l,u,v,d, &
                  txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                  qcxts5,qcyts6, &
                  equat,ncin, &
                  sn(npsn),lgsnlt, &
                  vol, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          else
!            schemas de Jameson standard et de Roe
             call sch_turb( &
                  l,u,v,d, &
                  txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
                  qcxts5,qcyts6, &
                  equat,ncin, &
                  sn(npsn),lgsnlt, &
                  vol, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz)
          endif
       elseif(kditur.eq.4) then
!           schema de Roe avec extrapolation MUSCL
!                call sch_turb_roe_muscl( &
!                 l,u,v, &
!                 txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
!                 qcxts5,qcyts6, &
!                 equat, &
!                 sn(npsn),lgsnlt, &
!                 vol,dvxx,dvxy,dvxz,dvyx)
       endif
!
!       ----------------------------------------------------------
!      Residu instationnaire DTS pour equation de transport
!
       if (kfmg.eq.3) then
          if(kdualns.eq.2) then
             call met_dual2( &
                  icycle,l,u,v, &
                  vol,ptdual)
          else
             call met_dual( &
                  icycle,l,u,v, &
                  vol,ptdual)
          endif
       endif
!------------------------------------------------------------------
!         phase implicite
!------------------------------------------------------------------
!
!         Jacobi par lignes alternes
       if(equat(3:4).eq.'2d') then
          call impli2_eqt( &
               l,u,dt,v, &
               mu,mut,cfke,ncin, &
               ncyc, &
               sn(npsn),lgsnlt, &
               vol, &
               dvxx,dvxy,dvxz,dvyx,dvyy, &
               dvyz,dvzx,dvzy,dvzz,dtdx)
       elseif(equat(3:4).eq.'3d') then
          call impli2_eqt_3d( &
               l,u,dt,v, &
               mu,mut,cfke,ncin, &
               ncyc, &
               sn(npsn),lgsnlt, &
               vol, &
               dvxx,dvxy,dvxz,dvyx,dvyy, &
               dvyz,dvzx,dvzy,dvzz,dtdx)
       endif
!
!          application des limiteurs
       call at_cutke(l,v)
!
    enddo  !fin de boucle sur les domaines
!
!*********************************************************************
!       calcul de la viscosite turbulente mut
!*********************************************************************
!
    do l=1,lzx
!
!------------------------------------------------------------------
!       Modele k-eps de Jones Launder ou Launder Sharma
!--------------------------------------------------------------------
!
       if(kcmut.eq.1) then
!
          if((equatt(4:4).eq.' ').or.(equatt(4:4).eq.'S').or. &
               (equatt(4:4).eq.'C').or.(equatt(4:4).eq.'L')) then
!            modele de base ou avec correction SST
             call met_kemut( &
                  l, &
                  sn,vol,t, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  dist,v,mu,mut, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
          else if(equatt(4:4).eq.'R') then
!            modele de Jones Launder realisable
             call met_kemutr( &
                  l,ncyc, &
                  sn,vol,t, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  v,mu,mut, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          else if(equatt(4:4).eq.'M') then
!            modele de Jones Launder modifie avec C_mu variable
             call met_kemutm( &
                  l, &
                  sn,vol,t, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  v,mu,mut, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          endif
!
!------------------------------------------------------------------
!       Modele de Spalart-Allmaras
!------------------------------------------------------------------
!
       else if(kcmut.eq.2) then
!
          call met_samut(l,v,mu,mut)
!
!------------------------------------------------------------------
!       Modele k-l de Smith
!------------------------------------------------------------------
!
       else if(kcmut.eq.3) then
!
          if((equatt(4:4).eq.' ').or.(equatt(4:4).eq.'L')) then
!            modele k-l de base et SAS
             call met_klmut(l,v,mu,mut,dist)
!
          elseif(equatt(4:4).eq.'S') then
!            modele k-l avec correction SST
             call met_klsmut( &
                  l, &
                  sn,vol,t, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  dist,v,mu,mut, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          elseif(equatt(4:4).eq.'R') then
!            modele k-l realisable
             call met_klrmut( &
                  l, &
                  sn,vol,t, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  dist,v,mu,mut, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
          else
             write(imp,'(/,''!!!met_num: modele k-l de Smith '',''non prevu'')')
             stop
          endif
!
!-------------------------------------------------------------
!       Modele k-w de Wilcox et de Menter
!-------------------------------------------------------------
!
       else if(kcmut.eq.4) then
!
          if(equatt(4:4).eq.'G') then
!         modele EASM de Gatski et Speziale
!            call met_komutg( &
!                 l, &
!                 sn,vol,t, &
!                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
!                 v,mu,mut, &
!                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          elseif((equatt(4:4).eq.'R').or. &
               (equatt(5:5).eq.'R')) then
!         Modeles k-omega de Menter avec realisabilite de Durbin
             call met_komutr( &
                  l, &
                  sn,vol,t, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  dist,v,mu,mut, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          else
!         Modeles k-omega de Wilcox et Menter
             call met_komut( &
                  l, &
                  sn,vol,t, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  dist,v,mu,mut, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          endif
!
!------------------------------------------------------------------
!       Modele k-eps de Chieng
!------------------------------------------------------------------
!
       else if(kcmut.eq.5) then
!
          call met_chmut( &
               l, &
               v,mu,mut,dist,mnpar,utau)
!
!------------------------------------------------------------------
!       Modele k-eps de Jones Launder ou Launder Sharma avec correction
!       pour utiliser la fonction f_mu de Smith dans les regions externes
!       et les sillages
!------------------------------------------------------------------
!
       else if(kcmut.eq.6) then
!
          call met_mutke2( &
               l,ncyc, &
               v,mu,mut,dist,mnpar,ncin, &
               txxf5x,txyf5y,txzf5z,tyyf6x,tyzf6y,tzzf6z, &
               qcz000)
!
!------------------------------------------------------------------
!       Modele bicouche k-eps RNG avec k-l a la paroi
!------------------------------------------------------------------
!
       else if(kcmut.eq.7) then
!
!            call met_rngmut( &
!                 l, &
!                 v,mu,mut,dist,fracmod)
!
!------------------------------------------------------------------
!       Modele k-omega de Kok et Wilcox compressible
!------------------------------------------------------------------
!
       else if(kcmut.eq.8) then
!
          if((equatt(4:4).eq.' ').or.(equatt(4:4).eq.'S')) then
!           modele standard ou SST
             call met_kocmut( &
                  l, &
                  sn,vol,t, &
                  dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
                  dist,v,mu,mut, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!            call met_kokmut(
!     &           l,
!     &           sn,vol,t,
!     &           dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz,
!     &           dist,v,mu,mut,
!     &           cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
          else if(equatt(4:4).eq.'R') then
!           modele avec realisabilite de Durbin
!            call met_kokmutr( &
!                 l, &
!                 sn,vol,t, &
!                 dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
!                 v,mu,mut, &
!                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
          else
             write(imp,'(/,''!!!met_num:  non prevu'')')
             stop
          endif
       endif
!
!------------------------------------------------------------------
!       Transition fixee
!------------------------------------------------------------------
!
       if(ktransi.eq.1) then
          call atctranske(l,v,mu,mut,mnpar,fgam)
       endif
    enddo  !fin de boucle sur les domaines

    DEALLOCATE(dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz,fracmod)
    deallocate(coefe)

!$OMP END MASTER
    return
  end subroutine met_num
end module mod_met_num
