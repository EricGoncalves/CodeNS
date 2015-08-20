module mod_smg_flu
  implicit none
contains
  subroutine smg_flu( &
       img,itypdf, &
       ff, &
       icyc,ncyc,icycle, &
       x,y,z,r,nxn,nyn,nzn, &
       sn, &
       vol, &
       tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
       u,mu,mut,dist, &
       v,d, &
       ptdual, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10, &
       tm11,tm12,tm13, &
       ncin,mnpar,fgam, &
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
!_DA  DATE_C : avril 2002 - Eric GONCALVES / SINUMEF
!
!     ACT
!_A    Calcul du bilan de flux et accumulation de "forcing function".
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use chainecarac
    use schemanum
    use maillage
    use boundary
    use modeleturb
    use proprieteflu
    use definition
    use mod_sch_ausmp
    use mod_sch_jameson3pond
    use mod_sch_rusanov
    use mod_sch_ausmp_prcd
    use mod_sch_weno5pond
    use mod_sch_jameson_pond
    use mod_sch_hllc_prcd
    use mod_dissip_jameson
    use mod_sch_dual
    use mod_sch_jameson3
    use mod_zvisqc
    use mod_sch_hllc_euler
    use mod_sch_weno3pond_3d
    use mod_sch_roe_pond
    use mod_sch_weno5_3d
    use mod_sch_jameson
    use mod_rbord
    use mod_sch_dual2
    use mod_sch_weno3pond
    use mod_sch_roe_prcd
    use mod_sch_weno3_3d
    use mod_rfve
    use mod_rfvc
    use mod_sch_roe_euler
    use mod_cllparoi1
    use mod_zpres
    use mod_sch_hllc
    use mod_sch_weno5
    use mod_sch_roe
    use mod_sch_weno3
    use mod_smg_res
    use mod_dissip_jameson_prcd2
    use mod_sch_ausmp_pond
    use mod_zvismo
    use mod_sch_rusanov_prcd
    implicit none
    integer          ::        icyc,     icycle,        img,     itypdf,          l
    integer          ::      lgsnlt,         lm,    mcyturb,        mfc,        mfn
    integer          ::   mnc(ip43),mnpar(ip12),  mnr(ip44), ncbd(ip41), ncin(ip41)
    integer          ::        ncyc,       npsn
    double precision ::  bceqt(ip41,neqt),      cmui1(ip21),      cmui2(ip21),      cmuj1(ip21),      cmuj2(ip21)
    double precision ::       cmuk1(ip21),      cmuk2(ip21),       cson(ip11),        cvi(ip21),        cvj(ip21)
    double precision ::         cvk(ip21),     d(ip11,ip60),        d0x(ip40),        d0y(ip40),        d0z(ip40)
    double precision ::        dist(ip12),    ff(ip11,ip60),       fgam(ip42),         mu(ip12),        mut(ip12)
    double precision ::         nxn(ip42),        nyn(ip42),        nzn(ip42),       pres(ip40),   pression(ip11)
    double precision :: ptdual(ip11,ip60),        qcx(ip12),        qcy(ip12),        qcz(ip12),          r(ip11)
    double precision ::         rod(ip40),       roed(ip40),       roud(ip40),       rovd(ip40),       rowd(ip40)
    double precision ::         rpi(ip40),        rti(ip40),    sn(ip31*ndir),        tm1(ip40),       tm10(ip40)
    double precision ::        tm11(ip40),       tm12(ip40),       tm13(ip40),        tm2(ip40),        tm3(ip40)
    double precision ::         tm4(ip40),        tm5(ip40),        tm6(ip40),        tm7(ip40),        tm8(ip40)
    double precision ::         tm9(ip40),        tn1(ip00),       tn10(ip00),        tn2(ip00),        tn3(ip00)
    double precision ::         tn4(ip00),        tn5(ip00),        tn6(ip00),        tn7(ip00),        tn8(ip00)
    double precision ::         tn9(ip00),       toxx(ip12),       toxy(ip12),       toxz(ip12),       toyy(ip12)
    double precision ::        toyz(ip12),       tozz(ip12),         tp(ip40),     u(ip11,ip60),     v(ip11,ip60)
    double precision ::         vol(ip11),          x(ip21),        xnr(ip44),          y(ip21),        ynr(ip44)
    double precision ::           z(ip21),        znr(ip44),      ztemp(ip11)
    double precision,allocatable :: dt(:)
    allocate(dt(ip11))
!
!-----------------------------------------------------------------------
!

!
    if (equat(1:2).eq.'ns') then
       if ( ncyturb .gt. 1 ) then
          mcyturb=mod(ncyc,ncyturb)
       else
          mcyturb = 1
       endif
    endif
!
    do l=1,lzx
       lm=l+(img-1)*lz
!
       call zpres( &
            lm,v, &
            pression,ztemp,cson)
!
!-------calcul de la viscosite moleculaire----------------------------
!
       if (equat(1:2).eq.'ns') then
          call zvismo(lm,mu,v,ztemp)
       endif
    enddo
!
!-------calcul de la dissipation artificielle Jameson----------------------
!
    if(ischema.le.4) then
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
       do l=1,lzx
          lm=l+(img-1)*lz
          npsn  =ndir*npfb(lm)+1
          lgsnlt=nnn(lm)
!
          if(kprec.eq.0) then
             call dissip_jameson( &
                  lm,v,d, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,pression,cson)

          elseif(kprec.eq.1) then
!               call dissip_jameson_prcd( &
!                 lm,v,d, &
!                 equat, &
!                 sn(npsn),lgsnlt, &
!                 tn2,tn3,pression,ztemp,cson)

          elseif(kprec.eq.2) then
             call dissip_jameson_prcd2( &
                  lm,v,d, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,pression,cson)
          endif
       enddo
    endif
!
!---remplissage des valeurs de v aux bords (centres des facettes frontieres)
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
!-----calcul des termes visqueux---------------------------------------------
!
    if (equat(1:2).eq.'ns') then
!
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
!
!-------application de la condition aux limites lois de paroi-------------------
!
       if((img.eq.mgl).and.(lparoi.eq.1)) then
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
       endif
    endif
!
!----application du schema numerique--------------------------------
!
    do l=1,lzx
!
       lm=l+(img-1)*lz
       npsn  =ndir*npfb(lm)+1
       lgsnlt=nnn(lm)
!
!********JAMESON**************************************************
       select case(ischema)
       case(1)
          call sch_jameson( &
               lm,0, &
               u,v,d,ff, &
               toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
               equat, &
               sn(npsn),lgsnlt, &
               tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
               pression)
!
       case(2)
          call sch_jameson_pond( &
               lm,0, &
               u,v,d,ff, &
               toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
               equat, &
               sn(npsn),lgsnlt, &
               tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
               pression, &
               cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
       case(3)
          call sch_jameson3( &
               lm,0, &
               u,v,d,ff, &
               toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
               equat, &
               sn(npsn),lgsnlt, &
               tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
               pression)
!
       case(4)
          call sch_jameson3pond( &
               lm,0, &
               u,v,d,ff, &
               toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
               equat, &
               sn(npsn),lgsnlt, &
               tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
               pression, &
               cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
               cvi,cvj,cvk)
!
!********AUSM+ de LIOU**************************************************
       case(5)
!
          if(equat(1:2).eq.'ns') then
             if(kprec.eq.0) then
                call sch_ausmp( &
                     lm,0, &
                     u,v,d,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)
             elseif(kprec.ge.1) then
                call sch_ausmp_prcd( &
                     lm,0, &
                     u,v,d,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)
             endif
          endif
!
!********AUSM+ pondere *************************************************
       case(6)
!
          if(equat(1:2).eq.'ns') then
             if(kprec.eq.0) then
                call sch_ausmp_pond( &
                     lm,0, &
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
!********ROE**************************************************
       case(7)
!
          if(equat(1:2).eq.'ns') then
             if(kprec.eq.0) then
                call sch_roe( &
                     lm,0, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)
!
             elseif(kprec.ge.1) then
                call sch_roe_prcd( &
                     lm,0, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)

             endif
          elseif(equat(1:2).eq.'eu') then
             call sch_roe_euler( &
                  lm,0, &
                  u,v,ff, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                  pression)
          endif
!
!--------Schema de Roe pondere------------------------------------
!
       case(8)
          if(kprec.eq.0) then
             call sch_roe_pond( &
                  lm,0, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                  pression, &
                  cvi,cvj,cvk, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
          elseif(kprec.ge.1) then
!           call sch_roe_pond_prcd( &
!                 lm,0, &
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
!--------Schema de Rusanov------------------------------------
!
       case(9)
          if(kprec.eq.0) then
             call sch_rusanov( &
                  lm,0, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                  pression)
          elseif(kprec.ge.1) then
             call sch_rusanov_prcd( &
                  lm,0, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                  pression)
          endif
!
!--------Schema de Jiang&Chu WENO ordre 3-----------------------------------
!
       case(10)
          if(equat(3:4).eq.'2d') then
             call sch_weno3( &
                  lm,0, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                  pression)
          elseif(equat(3:4).eq.'3d') then
             call sch_weno3_3d( &
                  lm,0, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                  pression)
          endif
!
!-------------- WENO 3 pondere-----------------------------------------
!
       case(11)
          if(equat(3:4).eq.'2d') then
             call sch_weno3pond( &
                  lm,0, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                  pression, &
                  cvi,cvj,cmui1,cmui2,cmuj1,cmuj2)
          elseif(equat(3:4).eq.'3d') then
             call sch_weno3pond_3d( &
                  lm,0, &
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
       case(12)
          if(equat(3:4).eq.'2d') then
             call sch_weno5( &
                  lm,0, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                  pression)
          elseif(equat(3:4).eq.'3d') then
             call sch_weno5_3d( &
                  lm,0, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                  pression)
          endif
!
!--------Schema WENO ordre 5 pondere----------------------------
!
       case(13)
          if(equat(3:4).eq.'2d') then
             call sch_weno5pond( &
                  lm,0, &
                  u,v,ff, &
                  toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                  pression, &
                  cvi,cvj,cmui1,cmui2,cmuj1,cmuj2)
          elseif(equat(3:4).eq.'3d') then
!              call sch_weno5pond_3d( &
!                 lm,0, &
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
!*********Schema HLLC avec extrapolation MUSCL*******************
!
       case(15)
          if(equat(1:2).eq.'ns') then
             if(kprec.eq.0) then
                call sch_hllc( &
                     lm,0, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)
             elseif(kprec.ge.1) then
                call sch_hllc_prcd( &
                     lm,0, &
                     u,v,ff, &
                     toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                     equat, &
                     sn(npsn),lgsnlt, &
                     tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                     pression)
             endif
          elseif(equat(1:2).eq.'eu') then
             call sch_hllc_euler( &
                  lm,0, &
                  u,v,ff, &
                  equat, &
                  sn(npsn),lgsnlt, &
                  tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                  pression)
          endif
       end select
!
!------calcul du residu instationnaire------------------------------------
!
       if(kfmg.eq.3) then
          if((kdualns.eq.0).or.(kdualns.eq.1)) then
             call sch_dual( &
                  lm,u,v,icycle, &
                  vol,ptdual)
          elseif(kdualns.eq.2) then
             call sch_dual2( &
                  lm,u,v,icycle, &
                  vol,ptdual)
          endif
       endif
!
!---calcul du bilan de flux et accumulation de "forcing function"-------
!
       call smg_res( &
            itypdf, &
            lm,u,ff)
!
    enddo
!
    deallocate(dt)
    return
  end subroutine smg_flu
end module mod_smg_flu
