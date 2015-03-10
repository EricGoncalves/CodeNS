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
!
!-----------------------------------------------------------------------
!
      real nxn,nyn,nzn
      real mu,mut
      dimension x(ip21),y(ip21),z(ip21)
      dimension u(ip11,ip60),v(ip11,ip60),d(ip11,ip60),ff(ip11,ip60)
      dimension ptdual(ip11,ip60)
      dimension sn(ip31*ndir)
      dimension vol(ip11),r(ip11),pression(ip11),ztemp(ip11),cson(ip11)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
      dimension bceqt(ip41,neqt)
      dimension rpi(ip40),rti(ip40)
      dimension d0x(ip40),d0y(ip40),d0z(ip40)
      dimension pres(ip40),tp(ip40)
      dimension rod(ip40),roud(ip40),rovd(ip40),rowd(ip40),roed(ip40)
      dimension xnr(ip44),ynr(ip44),znr(ip44),mnr(ip44)
      dimension mu(ip12),mut(ip12),toxx(ip12),toxy(ip12),toxz(ip12), &
                toyy(ip12),toyz(ip12),tozz(ip12),qcx(ip12),qcy(ip12), &
                qcz(ip12),dist(ip12)
      dimension tn1(ip00),tn2(ip00),tn3(ip00),tn4(ip00),tn5(ip00), &
                tn6(ip00),tn7(ip00),tn8(ip00),tn9(ip00),tn10(ip00)
      dimension tm1(ip40),tm2(ip40),tm3(ip40),tm4(ip40),tm5(ip40), &
                tm6(ip40),tm7(ip40),tm8(ip40),tm9(ip40), &
                tm10(ip40),tm11(ip40),tm12(ip40),tm13(ip40)
      dimension ncin(ip41),mnpar(ip12),fgam(ip42),mnc(ip43)
      dimension cvi(ip21),cvj(ip21),cvk(ip21),cmui1(ip21),cmui2(ip21), &
                cmuj1(ip21),cmuj2(ip21),cmuk1(ip21),cmuk2(ip21)
      dimension dt(ip11)

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
                pression,ztemp,cson)
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
         if(ischema.eq.1) then
          call sch_jameson( &
                 lm,0, &
                 u,v,d,ff, &
                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                 equat, &
                 sn(npsn),lgsnlt, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                 pression)
!
         elseif(ischema.eq.2) then
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
         elseif(ischema.eq.3) then
         call sch_jameson3( &
                 lm,0, &
                 u,v,d,ff, &
                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                 equat, &
                 sn(npsn),lgsnlt, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9, &
                 pression)
!
         elseif(ischema.eq.4) then
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
         elseif(ischema.eq.5) then
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
         elseif(ischema.eq.6) then
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
         elseif(ischema.eq.7) then
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
         elseif(ischema.eq.8) then
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
         elseif(ischema.eq.9) then
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
         elseif(ischema.eq.10) then
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
         elseif(ischema.eq.11) then
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
         elseif(ischema.eq.12) then
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
         elseif(ischema.eq.13) then
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
         elseif(ischema.eq.15) then
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
        endif
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
      return
      end
