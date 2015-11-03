module mod_atsch_num
  implicit none
  integer         ,allocatable ::  m1tb(:), m2tb(:),nfrtb(:)
contains
  subroutine atsch_num( &
       img,ityprk, &
       icyc,ncyc,idcyc,icycle, &
       mg,  &
       x,y,z,r,nxn,nyn,nzn, &
       sn, &
       vol, &
       tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
       u,mu,mut,dist,cfke, &
       mnpar,fgam,utau, &
       v,dt,d, &
       ptdual, &
       ff, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10, &
       tm11,tm12,tm13, &
       ncin, &
       mnc, &
       ncbd,mnr,xnr,ynr,znr, &
       bceqt, &
       rpi,rti,d0x,d0y,d0z,pres,tp, &
       rod,roud,rovd,rowd,roed, &
       pression,ztemp,cson, &
       cvi,cvj,cvk, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DA  DATE_C : mars 2002 - Eric GONCALVES / SINUMEF
!
!     ACT
!_A    Realisation d'une iteration du schema numerique sur un
!_A    niveau de grille donne pour toute la configuration de calcul.
!
!     INP
!_I    icyc       : arg int    ; cycle courant du calcul (initialise
!_I                                        a numt dans "flec")
!_I    ncyc       : arg int    ; cycle courant de l'execution (remis a
!_I                                        a 1 a chaque reprise)
!_I    keinit     : com int    ; controle initialisation et calcul
!_I                                        equations de transport
!_I    icytur0    : arg int    ; nbr de cycle en deb de calcul au cours
!_I                              desquelles mut n'est pas mis a jour
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use chainecarac
    use sortiefichier
    use schemanum
    use modeleturb
    use maillage
    use boundary
    use definition
    use proprieteflu
    use mod_atcaldis
    use mod_atcaldis3
    use mod_atintrans
    use mod_cllparoi1
    use mod_cllparoi2
    use mod_chrono
    use mod_dissip_jameson
    use mod_dissip_jameson_turb
    use mod_dissip_jameson_prcd2
    use mod_sch_dual
    use mod_sch_dual2
    use mod_sch_expli
    use mod_implimf_2d
    use mod_implimf_3d
    use mod_implimf_prcd2
    use mod_implimf_eu
    use mod_met_num
    use mod_met_ini
    use mod_met_inmut
    use mod_met_inisa
    use mod_met_inikl
    use mod_met_iniko
    use mod_met_iniuttau
    use mod_met_rfve
    use mod_met_rfvc
    use mod_prcd_turkel
    use mod_rbord
    use mod_rfvc
    use mod_rfve
    use mod_rbtc
    use mod_sch_acou
    use mod_sch_ausmp
    use mod_sch_ausmp_prcd
    use mod_sch_ausmp_pond
    use mod_sch_jameson
    use mod_sch_jameson3
    use mod_sch_jameson_pond
    use mod_sch_jameson3pond
    use mod_sch_hllc
    use mod_sch_hllc_prcd
    use mod_sch_hllc_euler
    use mod_sch_roe
    use mod_sch_roe_euler
    use mod_sch_roe_pond
    use mod_sch_roe_prcd
    use mod_sch_rusanov
    use mod_sch_rusanov_prcd
    use mod_sch_weno3
    use mod_sch_weno3split
    use mod_sch_weno3pond_split2
    use mod_sch_weno3pond_3d
    use mod_sch_weno3pond
    use mod_sch_weno3split2
    use mod_sch_weno3_3d
    use mod_sch_weno5pond
    use mod_sch_weno5_3d
    use mod_sch_weno5
    use mod_zpres
    use mod_zvisqc
    use mod_mpi
    implicit none
    integer          ::        icyc,     icycle,      idcyc,        img,     ityprk
    integer          ::           l,     ldismx,     lgsnlt,         lm,    mcychro
    integer          ::     mcyturb,        mfc,        mfn,         mg,  mnc(ip43)
    integer          :: mnpar(ip12),  mnr(ip44),         nc, ncbd(ip41), ncin(ip41)
    integer          ::        ncyc,       ndeb,       nfin,       npsn
    double precision ::  bceqt(ip41,neqt),       cfke(ip13),      cmui1(ip21),      cmui2(ip21),      cmuj1(ip21)
    double precision ::       cmuj2(ip21),      cmuk1(ip21),      cmuk2(ip21),       cson(ip11),        cvi(ip21)
    double precision ::         cvj(ip21),        cvk(ip21),     d(ip11,ip60),        d0x(ip40),        d0y(ip40)
    double precision ::         d0z(ip40),       dist(ip12),         dt(ip11),            dtpas,    ff(ip11,ip60)
    double precision ::        fgam(ip42),         mu(ip12),        mut(ip12),        nxn(ip42),        nyn(ip42)
    double precision ::         nzn(ip42),       pres(ip40),   pression(ip11),ptdual(ip11,ip60),        qcx(ip12)
    double precision ::         qcy(ip12),        qcz(ip12),          r(ip11),        rod(ip40),       roed(ip40)
    double precision ::        roud(ip40),       rovd(ip40),       rowd(ip40),        rpi(ip40),        rti(ip40)
    double precision ::     sn(ip31*ndir),        tm1(ip40),       tm10(ip40),       tm11(ip40),       tm12(ip40)
    double precision ::        tm13(ip40),        tm2(ip40),        tm3(ip40),        tm4(ip40),        tm5(ip40)
    double precision ::         tm6(ip40),        tm7(ip40),        tm8(ip40),        tm9(ip40),        tn1(ip00)
    double precision ::        tn10(ip00),        tn2(ip00),        tn3(ip00),        tn4(ip00),        tn5(ip00)
    double precision ::         tn6(ip00),        tn7(ip00),        tn8(ip00),        tn9(ip00),       toxx(ip12)
    double precision ::        toxy(ip12),       toxz(ip12),       toyy(ip12),       toyz(ip12),       tozz(ip12)
    double precision ::          tp(ip40),     u(ip11,ip60),       utau(ip42),     v(ip11,ip60),        vol(ip11)
    double precision ::           x(ip21),        xnr(ip44),          y(ip21),        ynr(ip44),          z(ip21)
    double precision ::         znr(ip44),      ztemp(ip11)
    logical          :: gfetke

    allocate(m2tb(ip00))
    allocate(nfrtb(ip00))
    allocate(m1tb(ip00))
!
!-----------------------------------------------------------------------
!
!


!     keinit : controle initialisation variables k et seconde (epsilon, l, ...)
!              avant "atsch_num" :
!                0      : lecture mut, k, epsilon
!                1 ou 2 : lecture mut, k et epsilon mis a 0
!              dans "atsch_num" :
!                0      : integration des equations de transport
!                1      : initialisation k et epsilon a partir du mut lu
!                2      : initialisation k et epsilon a partir de mut=mu
!              keinit passe de 1 ou 2 a 0 en meme temps que k et epsilon
!              sont initialise ce qui se produit pour
!              "icyc" sup ou egal "icytur0" et "keinit" non nul
!
!     grille fine et calcul equations de transport : gfetke=.true.
!      gfetke=(equat(6:7).eq.'ke') .and. (img.eq.mg)
    gfetke=(equat(6:7).eq.'ke') .and. (img.eq.1)
!
    dtpas = dt1min
    if(ncychro.gt.1) then
       mcychro=mod(ncyc,ncychro)
    else
       mcychro = 1
    endif
!
    if (equat(1:2).eq.'ns') then
       if ( ncyturb .gt. 1 ) then
!        ncyturb: frequence integration equations de transport (fatdon)
          mcyturb=mod(ncyc,ncyturb)
       else
          mcyturb = 1
       endif
    endif
!
!-----------------------------------------------------------
!     calcul de la distance aux parois
!-----------------------------------------------------------
!
    if (ncyc.eq.1 .and. img.eq.mg ) then
!        calcul distance pour epaisseurs integrales si necessaire
       SELECT CASE(kcaldis)
       CASE(1:2)
!          1->integration suivant lignes de maillage
!          2->minimum distance cellule-facette sans optimisation
          call atcaldis( &
               x,y,z,nxn,nyn,nzn, &
               tn1,tn2,tn3,tn4,tn5,tn6,tn7, &
               dist,mnpar,fgam,img, &
               ncin,mnc,ncbd)
!
       CASE(3)
!          3->minimum distance cellule-facette avec optimisation
          call atcaldis3( &
               x,y,z,nxn,nyn,nzn, &
               tn1,tn2,tn3,tn4,tn5,tn6,tn7, &
               dist,mnpar,fgam, &
               ncin,mnc,ncbd, &
               m1tb,m2tb,nfrtb)
!
       CASE(0)
          if(klecdis.eq.1) then
!          lecture des distances
             call at_lecdist( &
                  ldismx,dist,mnpar)

          endif
       case DEFAULT
          write(imp,'(/,''!atsch_num calcul distance non prevu STOP'')')
          stop
       END SELECT
       if(kecrdis.eq.1) then
!          ecriture "fdist" des distances pour tous les domaines
          call at_ecrdist( &
               0, &
               dist,mnpar)
       endif
!
!--------initialisation fonction "fgam" pour transition. Lecture "fatdon"
!
       call atintrans(ncin,fgam)
    endif
!
!-------------------------------------------------------------------
!     calcul des variables thermo
!-------------------------------------------------------------------
!
    do l=1,lzx
       lm=l+(img-1)*lz
       call zpres( &
            lm,v, &
            pression,ztemp,cson)
    enddo
!
!-----------------------------------------------------------
!     initialisation de k et epsilon
!-----------------------------------------------------------
!
    if(gfetke) then
       if(kinke.eq.3) then
          alfak=sigmak
          alfae=sigmal
       endif
!
       if(keinit.ge.1) then
          do l=1,lzx
             lm = l+(img-1)*lz
             call zvismo(lm,mu,v,ztemp)
!
             ndeb=npc(lm)+1
             nfin=npc(lm)+nnc(lm)
             do nc=ndeb,nfin
                v(nc,6)=epsk
                v(nc,7)=epse
             enddo
          enddo
!
          if(keinit.eq.2) then
!        mut mis egal a mu/10
             do l=1,lzx
                call met_inmut(l,mu,mut)
             enddo
          endif
!
          if(icyc.ge.icytur0) then
!
             do mfn=1,mtnx
                lbd(mfn)=nfbn(mfn)+(img-1)*mtb
             enddo
             nbd=mtnx
             call rfve( &
                  v,pression,ztemp,cson, &
                  ncbd,ncin)
!
             do mfc=1,mtcx
                lbd(mfc)=nfbc(mfc)+(img-1)*mtb
             enddo
             nbd=mtcx
             call rfvc( &
                  v,ncbd,mnc, &
                  pression,ztemp,cson,ncin)
!
              do mfc=1,mtcx
                 lbd(mfc)=nfbc(mfc)
              enddo
              nbd=mtcx
              call rbtc( & ! TODO: before or after ?
                   toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                   ncbd,ncin,mnc)

             do l=1,lzx
                SELECT CASE(kinke)
                CASE(1)
!
                   call met_ini( &
                        l,v,mut,mu, &
                        sn,vol,tn1, &
                        tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                        cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
                   if(kutau.eq.1) then
!
!             calcul du tenseur des contraintes visqueuses pour
!             le calcul de "utau" sur toutes les parois (debut de "met_num")
!             modeles k-omega et modele de Chien
                      call met_iniuttau( &
                           l,mu,mut,v,equat, &
                           toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                           sn,vol, &
                           ncbd,ncin,mnc, &
                           mnr,xnr,ynr,znr, &
                           nxn,nyn,nzn, &
                           tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                           cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
                           ztemp)
!
                   endif
!
                CASE(2)
!
!------------Modele de Spalart Allmaras------------------------------
!
                   call met_inisa(l,v,mut,mu)
!
                CASE(3)
!
!-------------Modele k-l de Smith----------------------------------
!
                   call met_inikl( &
                        l,v,mut,mu,dist, &
                        sn,vol,tn1, &
                        tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                        cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
                CASE(4)
!
!-------------Modele k-omega de Wilcox ou Menter----------------------
!             Modele k-omega bas Reynolds de Wilcox
!             Modele k-omega de Wilcox ou Menter avec rugosite
!
                   if(kutau.eq.1) then
!
!             calcul du tenseur des contraintes visqueuses pour
!             le calcul de "utau" sur toutes les parois (debut de "met_num")
!             modeles k-omega et modele de Chien
!
                      call met_iniuttau( &
                           l,mu,mut,v,equat, &
                           toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                           sn,vol, &
                           ncbd,ncin,mnc, &
                           mnr,xnr,ynr,znr, &
                           nxn,nyn,nzn, &
                           tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                           cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
                           ztemp)
!
                   endif
!
                   call met_iniko( &
                        l,ncin,ncbd, &
                        v,mut,mu,dist,mnpar, &
                        sn,vol,tn1, &
                        tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                        cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
                CASE(6)
!
                   write(imp,'(/,"!!!atsch_num: kinke=",i4,4x,"commencer le calcul avec k-l de Smith")')kinke
                   stop
!
                CASE DEFAULT
                   write(imp,'(/,"!!!atsch_num: kinke=",i4,4x,"non prevue-STOP")')kinke
                   stop
                END SELECT
             enddo   ! fin boucle domaine
!
              do mfc=1,mtcx
                 lbd(mfc)=nfbc(mfc)
              enddo
              nbd=mtcx
              call rbtc( & ! TODO: before or after ?
                   toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                   ncbd,ncin,mnc)

             if(kutau.eq.1) then
!
!           initialisation de "utau" pour modele de Chien, le
!           modele k-omega bas Reynolds de Wilcox et le modele
!           k-omega de Wilcox ou Menter avec rugosite
                call met_uttau( &
                     ncbd, &
                     nxn,nyn,nzn, &
                     toxx,toxy,toxz,toyy,toyz,tozz, &
                     v,utau)
!
             endif
             if((kfmg.ne.2).or.(img.eq.1))  keinit=0
          endif      !fin test sur icytur0
       endif        !fin test sur keinit superieur ou egal a 1
    endif          !fin test sur gfetke
!
!******************************************************************************
!    calcul du point courant par le schema numerique (tous niveaux de grilles)
!******************************************************************************
!
    do l=1,lzx
       lm=l+(img-1)*lz
!
!-------calcul de la viscosite moleculaire------------------------------------
!
       if (equat(1:2).eq.'ns') then
          call zvismo(lm,mu,v,ztemp)
       endif
    enddo
!
!-------calcul du pas de temps------------------------------------------------
!
    if((mcychro.eq.1).or.(icyc.lt.icychr0)) then
       call chrono( &
            img, &
            mu,mut,v,dt,dtpas, &
            sn,vol, &
            tn1,tn2,tn3,tn4, &
            cson,pression)
    endif
!
!-------calcul dissipation artificielle du schema de Jameson--------------
!
!      champ moyen
    if(ischema.le.4) then

!-------prolongement des variables aux bords-------------------------------
    do mfn=1,mtnx
      lbd(mfn)=nfbn(mfn)+(img-1)*mtb
    enddo
    nbd=mtnx
    if(ncyc.eq.1) then
      call rfve( &
                 v,pression,ztemp,cson, &
                 ncbd,ncin)
    endif
    do mfc=1,mtcx
     lbd(mfc)=nfbc(mfc)+(img-1)*mtb
    enddo
    nbd=mtcx
    call rfvc( &
         v,ncbd,mnc, &
         pression,ztemp,cson,ncin)

       do l=1,lzx
          lm=l+(img-1)*lz
          npsn  =ndir*npfb(lm)+1
          lgsnlt=nnn(lm)
!
          SELECT CASE(kprec)
          CASE(0)
             call dissip_jameson( &
                  lm,v,d, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,pression,cson)
          CASE(1)  !(P,u,S)
!            call dissip_jameson_prcd( &
!                 lm,v,d, &
!                 equat, &
!                 sn(npsn),lgsnlt, &
!                 tn2,tn3,pression,ztemp,cson)
          CASE(2)  !(P,u,e)
             call dissip_jameson_prcd2( &
                  lm,v,d, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,pression,cson)
          END SELECT
       enddo
    endif
!
!       champ turbulent
    if((kditur.eq.1).or.(kditur.eq.3)) then
       do mfn=1,mtnx
          lbd(mfn)=nfbn(mfn)+(img-1)*mtb
       enddo
       nbd=mtnx
       call met_rfve(v,ncbd,ncin)
!
       do mfc=1,mtcx
          lbd(mfc)=nfbc(mfc)+(img-1)*mtb
       enddo
       nbd=mtcx
       call met_rfvc(v,ncbd,mnc,ncin)
!
       do l=1,lzx
          lm=l+(img-1)*lz
          npsn  =ndir*npfb(lm)+1
          lgsnlt=nnn(lm)
          call dissip_jameson_turb( &
               lm,v,d, &
               equat, &
               sn(npsn),lgsnlt, &
               tn1,pression)
       enddo
    endif
!
!-------------------------------------------------------------------------------
!         Continuite du champ moyen et application des conditions limites
!         pour calcul du flux verifiant les conditions limites
!-------------------------------------------------------------------------------
!
    call rbord( &
         0,img, &
         mnc,ncin,mnr,xnr,ynr,znr,nxn,nyn,nzn,ncbd, &
         sn,vol,v,dt, &
         bceqt, &
         rpi,rti,d0x,d0y,d0z,x,y,z, &
         pres,tp,rod,roud,rovd,rowd,roed, &
         icyc, &
         tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10,tm11, &
         tm12,tm13, &
         pression,ztemp,cson)
!
!--------------------------------------------------------------------------
!         integration du systeme turbulent
!--------------------------------------------------------------------------
!
    if (gfetke .and. keinit.eq.0 .and. ncyc.ge.icytur0) then
!
       call met_num( &
            ncbd,ncin,mnc,ncyc, &
            mnr,xnr,ynr,znr, &
            bceqt, &
            u,v,dt,d,mut,mu,cfke, &
            mnpar,dist,fgam,utau,r, &  !r utilise pour topz
            nxn,nyn,nzn, &
            toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
            sn,vol, &
            icycle,ptdual, &
            x,y,z, &
            tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
            tp,cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
            pression,ztemp,cson)
    endif
!
!-------calcul du tenseur des contraintes et du flux de chaleur------ ------
!
    if (equat(1:2).eq.'ns') then
       call zvisqc( &
            img, &
            v,mu,v(1,1), &
            x,y,z,mut,dist,mnpar,fgam, &
            toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
            icyc,mcyturb, &
            ncbd,ncin,mnc, &
            sn,vol, &
            tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
            cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
            ztemp)
    endif
!
!-------application de la condition aux limites lois de paroi----------------
!
    if(lparoi.eq.1) then
!        test pour geometrie avec coin!
!        if((lparoi.eq.1).and.(ncyc.ge.icytur0)) then
       if(img.eq.mg) then
          call cllparoi1( &
               img,ncyc, &
               v,dist, &
               mu,mut, &
               nxn,nyn,nzn, &
               ncin,ncbd, &
               toxx,toxy,toxz,toyy,toyz,tozz, &
               qcx,qcy,qcz, &
               mnpar,fgam,tp, &
               ztemp)
       else
          call cllparoi1( &
               1,ncyc, &
               v,dist, &
               mu,mut, &
               nxn,nyn,nzn, &
               ncin,ncbd, &
               toxx,toxy,toxz,toyy,toyz,tozz, &
               qcx,qcy,qcz, &
               mnpar,fgam,tp, &
               ztemp)
       endif
    elseif(lparoi.eq.2) then
!         approche de Smith
       call cllparoi2( &
            img,ncyc, &
            v, &
            nxn,nyn,nzn, &
            ncin,ncbd, &
            toxx,toxy,toxz,toyy,toyz,tozz, &
            qcx,qcy,qcz, &
            ztemp,utau,r)
    endif
!
!----------------------------------------------------------------------
!       calcul de l'increment explicite (calcul des flux numeriques)
!----------------------------------------------------------------------
!
    do l=1,lzx
       lm=l+(img-1)*lz
       npsn  =ndir*npfb(lm)+1
       lgsnlt=nnn(lm)
!
!--------Schema de Jameson--------------------------------------------
!
       SELECT CASE(ischema)
       CASE(1)
          call sch_jameson( &
               lm,ityprk, &
               u,v,d,ff, &
               toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
               equat, &
               sn(npsn),lgsnlt, &
               tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
               pression)
!
       CASE(2)
!         version ponderee
          call sch_jameson_pond( &
               lm,ityprk, &
               u,v,d,ff, &
               toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
               equat, &
               sn(npsn),lgsnlt, &
               tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
               pression, &
               cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
       CASE(3)
!         version ordre 3 avec correction de l'erreur dispersive
          call sch_jameson3( &
               lm,ityprk, &
               u,v,d,ff, &
               toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
               equat, &
               sn(npsn),lgsnlt, &
               tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
               pression)
!
       CASE(4)
!         version ordre 3 avec correction de l'erreur dispersive + ponderation
          call sch_jameson3pond( &
               lm,ityprk, &
               u,v,d,ff, &
               toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
               equat, &
               sn(npsn),lgsnlt, &
               tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
               pression, &
               cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
               cvi,cvj,cvk)
!
!--------Schema AUSM+ de Liou--------------------------------------------
!
       CASE(5)
          if(equat(1:2).eq.'ns') then
             if(kprec.eq.0) then
                call sch_ausmp( &
                     lm,ityprk, &
                     u,v,d,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)
             elseif(kprec.ge.1) then
                call sch_ausmp_prcd( &
                     lm,ityprk, &
                     u,v,d,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)
             endif
          endif
!
!--------Schema AUSM+ pondere--------------------------------------
!
       CASE(6)
          if(equat(1:2).eq.'ns') then
             if(kprec.eq.0) then
                call sch_ausmp_pond( &
                     lm,ityprk, &
                     u,v,d,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression, &
                     cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
                     cvi,cvj,cvk)
             endif
          endif
!
!--------Schema de Roe--------------------------------------------
!
       CASE(7)
!
          if(equat(1:2).eq.'ns') then
             if(kprec.eq.0) then
                call sch_roe( &
                     lm,ityprk, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)
!
             elseif(kprec.ge.1) then
                call sch_roe_prcd( &
                     lm,ityprk, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)
             endif
          elseif(equat(1:2).eq.'eu') then
             call sch_roe_euler( &
                  lm,ityprk, &
                  u,v,ff, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                  pression)
          endif
!
!--------Schema de Roe pondere------------------------------------
!
       CASE(8)
          if(kprec.eq.0) then
             call sch_roe_pond( &
                  lm,ityprk, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                  pression, &
                  cvi,cvj,cvk, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          elseif(kprec.ge.1) then
!             call sch_roe_pond_prcd( &
!                 lm,ityprk, &
!                 u,v,ff, &
!                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
!                 equat, &
!                 sn(npsn),lgsnlt, &
!                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
!                 pression, &
!                 cvi,cvj,cvk, &
!                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          endif
!
!--------Schema de Rusanov avec extrapolation MUSCL----------------------------------
!
       CASE(9)
!
          if(kprec.eq.0) then
             call sch_rusanov( &
                  lm,ityprk, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                  pression)
          elseif(kprec.ge.1) then
             call sch_rusanov_prcd( &
                  lm,ityprk, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                  pression)
          endif
!--------Schema de Jiang&Chu WENO ordre 3-----------------------------------
!
       CASE(10)
!           if(kprec.eq.0) then
          if(equat(3:4).eq.'2d') then
! attention: muscl sert de cle pour splitting
             SELECT CASE(muscl)
             CASE(0)
                call sch_weno3( &
                     lm,ityprk, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                     pression)
             CASE(1) !Steger&Warming
                call sch_weno3split( &
                     lm,ityprk, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                     pression)
             CASE(2) !Lax&Friedrichs
                call sch_weno3split2( &
                     lm,ityprk, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                     pression)
             END SELECT
          elseif(equat(3:4).eq.'3d') then
             call sch_weno3_3d( &
                  lm,ityprk, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                  pression)
!           elseif(kprec.ge.1) then
!                 call sch_weno3_prcd( &
!                 lm,ityprk, &
!                 u,v,ff, &
!                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
!                 equat, &
!                 sn(npsn),lgsnlt,
!                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
!                 pression,cson)
          endif
!
!--------Schema WENO ordre 3 pondere-----------------------------------
!
       CASE(11)
          if(equat(3:4).eq.'2d') then
! attention: muscl sert de cle pour splitting
             if(muscl.eq.0) then
                call sch_weno3pond( &
                     lm,ityprk, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                     pression, &
                     cvi,cvj,cmui1,cmui2,cmuj1,cmuj2)
             elseif(muscl.eq.2) then
                call sch_weno3pond_split2( &
                     lm,ityprk, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                     pression, &
                     cvi,cvj,cmui1,cmui2,cmuj1,cmuj2)
             endif
          elseif(equat(3:4).eq.'3d') then
             call sch_weno3pond_3d( &
                  lm,ityprk, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                  pression, &
                  cvi,cvj,cvk, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          endif
!
!--------Schema WENO ordre 5-----------------------------------
!
       CASE(12)
          if(equat(3:4).eq.'2d') then
             call sch_weno5( &
                  lm,ityprk, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                  pression)
          elseif(equat(3:4).eq.'3d') then
             call sch_weno5_3d( &
                  lm,ityprk, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                  pression)
          endif
!
!--------Schema WENO ordre 5 pondere---------------------------------
!
       CASE(13)
          if(equat(3:4).eq.'2d') then
             call sch_weno5pond( &
                  lm,ityprk, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                  pression, &
                  cvi,cvj,cmui1,cmui2,cmuj1,cmuj2)
          elseif(equat(3:4).eq.'3d') then
!              call sch_weno5_pond_3d( &
!                 lm,ityprk, &
!                 u,v,ff, &
!                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
!                 equat, &
!                 sn(npsn),lgsnlt, &
!                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
!                 pression, &
!                 cvi,cvj,cvk, &
!                cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          endif
!
!--------Schema HLLC avec extrapolation MUSCL-------------------------
!
       CASE(14)
!
          if(equat(1:2).eq.'ns') then
             if(kprec.eq.0) then
                call sch_hllc( &
                     lm,ityprk, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)
             elseif(kprec.ge.1) then
                call sch_hllc_prcd( &
                     lm,ityprk, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)
             endif
          else
             call sch_hllc_euler( &
                  lm,ityprk, &
                  u,v,ff, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                  pression)
          endif
       END SELECT
!
!--------Preconditionnement basse vitesse de Turkel--------------------
!        calcul des residus preconditionnes pour calcul tout explicite
!
       if((kprec.ge.1).and.(kmf(lm).eq.0)) then
          call prcd_turkel( &
               lm,u,v, &
               sn(npsn),lgsnlt, &
               ztemp,cson)
       endif
!
!--------calcul du residu instationnaire------------------------------
!
       if(kfmg.eq.3) then
          if((kdualns.eq.0).or.(kdualns.eq.1)) then !ordre 2 en temps
             call sch_dual( &
                  lm,u,v,icycle, &
                  vol,ptdual)
          elseif(kdualns.eq.2) then   !ordre 3 en temps
             call sch_dual2( &
                  lm,u,v,icycle, &
                  vol,ptdual)
          endif
       endif
!
!--------contribution acoustique ----------------------------------------
!
       if(lacou.eq.1) then
          call sch_acou( &
               lm,ityprk, &
               u,ff,dtpas, &
               x,y,z,idcyc, &
               vol)
       endif
!
!*************************************************************************
!      avance d'un pas de temps des variables
!*************************************************************************
!
       if(kmf(lm).eq.0) then
!        calcul tout explicite
          call sch_expli(lm,u,v,dt,vol)
!
       elseif(kmf(lm).eq.1) then
!-------------------------------------------------------------------------
!        phase implicite Matrix Free - Jacobi par points
!-------------------------------------------------------------------------
!
          if(equat(1:2).eq.'ns') then
             if(kprec.eq.0) then
                if(equat(3:4).eq.'2d') then
                   call implimf_2d( &
                        lm,u,dt,v,d,ff, &
                        mu,mut, &
                        equat,lmax(lm), &
                        sn(npsn),lgsnlt, &
                        vol,dtpas,ityprk, &
                        tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                        pression,cson)
                elseif(equat(3:4).eq.'3d') then
                   call implimf_3d( &
                        lm,u,dt,v,d,ff, &
                        mu,mut, &
                        equat,lmax(lm), &
                        sn(npsn),lgsnlt, &
                        vol,dtpas,ityprk, &
                        tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                        pression,cson)
                endif
!          elseif(kprec.eq.1) then
!         preconditionnement basse vitesse de Turkel (P,u,S)
!                call implimf_prcd( &
!                 lm,u,dt,v,d,ff, &
!                 mu,mut, &
!                 equat,lmax(lm), &
!                 sn(npsn),lgsnlt, &
!                 vol,dtpas,ityprk, &
!                 pression,ztemp,cson)
             elseif(kprec.eq.2) then
!         preconditionnement basse vitesse de Turkel (P,u,e)
!            if(equat(3:4).eq.'2d') then
                call implimf_prcd2( &
                     lm,u,dt,v,d,ff, &
                     mu,mut, &
                     equat,lmax(lm), &
                     sn(npsn),lgsnlt, &
                     vol,dtpas,ityprk, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression,cson)
!             elseif(equat(3:4).eq.'3d') then
!              call implimf_prcd2_3d( &
!                 lm,u,dt,v,d,ff, &
!                 mu,mut, &
!                 equat,lmax(lm), &
!                 sn(npsn),lgsnlt, &
!                 vol,dtpas,ityprk, &
!                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
!                 pression,cson)
!             endif
             endif
!             if(kdtl.eq.0) then
!              call implimf_prcd2inst( &
!                 lm,u,dt,v,d,ff, &
!                 mu,mut, &
!                 equat,lmax(lm), &
!                 sn(npsn),lgsnlt, &
!                 vol,dtpas,ityprk, &
!                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
!                 pression,cson)
!             endif
!
!-------Equations d'Euler------------------------------------------
!
          elseif(equat(1:2).eq.'eu') then
!          if(kprec.eq.0) then
             call implimf_eu( &
                  lm,u,dt,v,d,ff, &
                  equat,lmax(lm), &
                  sn(npsn),lgsnlt, &
                  vol,dtpas,ityprk, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                  pression,cson)
!
!          elseif(kprec.eq.1) then
!         preconditionnement basse vitesse de Turkel (P,u,s)
!            call implimf_eu_prcd( &
!                 lm,u,dt,v,d,ff, &
!                 equat,lmax(lm), &
!                 sn(npsn),lgsnlt, &
!                 vol,dtpas,ityprk, &
!                 pression,ztemp,cson)
!          elseif(kprec.eq.2) then
!         preconditionnement basse vitesse de Turkel (P,u,e)
!              call implimf_eu_prcd2( &
!                 lm,u,dt,v,d,ff, &
!                 equat,lmax(lm), &
!                 sn(npsn),lgsnlt, &
!                 vol,dtpas,ityprk, &
!                 pression,cson)
!          endif
          endif
!
!-------------------------------------------------------------------------
!        phase implicite bloc - Jacobi par lignes
!-------------------------------------------------------------------------
       elseif(kmf(lm).eq.2) then
          if(equat(1:2).eq.'eu') then
!           call impli_bloc_euler( &
!                 lm,u,dt,v, &
!                 equat,lmax(lm), &
!                 sn(npsn),lgsnlt, &
!                 vol,ityprk, &
!                 pression,ztemp,cson)
          endif
       endif
!
    enddo      !fin boucle sur domaines

    deallocate(m2tb,nfrtb,m1tb)
    return
  end subroutine atsch_num
end module mod_atsch_num
