module mod_smg_num
use mod_smg_cn
use mod_rbc
use mod_smg_cf
implicit none
contains
      subroutine smg_num( &
                 mg,ncycl, &
                 idcyc, &
                 x,y,z,r,nxn,nyn,nzn, &
                 sn, &
                 vol, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                 mu,mut,dist,cfke, &
                 mnpar,fgam,utau, &
                 v,dt, &
                 ptdual,vdual,vdual1,vdual2, &
                 d_volt,u_duv,ff_du,u0, &
                 toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
                 tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10, &
                 tm11,tm12,tm13, &
                 ncin, &
                 mnc, &
                 ncbd,mnr,xnr,ynr,znr, &
                 bceqt, &
                 rpi,rti,d0x,d0y,d0z,qtx,qty,qtz,pres,tp, &
                 rod,roud,rovd,rowd,roed, &
                 pression,ztemp,cson, &
                 cvi,cvj,cvk, &
                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!_DA  DATE_C : novembre 2001 - Eric Goncalves / SINUMEF
!
!     ACT
!_A    Realisation de ncycl cycles multigrille sur un niveau
!_A    donne de maillage emboite.
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use proprieteflu
      use schemanum
      use boundary
      use sortiefichier
      use definition
      use modeleturb
      use chainecarac
use mod_smg_fcr
use mod_sch_duup2
use mod_dua_resro
use mod_smg_upc
use mod_sortieplot
use mod_smg_flu
use mod_utitfr
use mod_sortieplot2
use mod_sch_duin
use mod_svfw
use mod_sortietest
use mod_smg_fcv
use mod_sch_duup
use mod_rbc
use mod_rscpsv
use mod_atsch_num
implicit none
integer :: mg
double precision :: x
double precision :: y
double precision :: z
double precision :: r
double precision :: sn
double precision :: vol
double precision :: tn1
double precision :: tn2
double precision :: tn3
double precision :: tn4
double precision :: tn5
double precision :: tn6
double precision :: tn7
double precision :: tn8
double precision :: tn9
double precision :: tn10
double precision :: dist
double precision :: cfke
integer :: mnpar
double precision :: fgam
double precision :: utau
double precision :: v
double precision :: dt
double precision :: ptdual
double precision :: vdual
double precision :: vdual1
double precision :: vdual2
double precision :: d_volt
double precision :: u_duv
double precision :: ff_du
double precision :: u0
double precision :: toxx
double precision :: toxy
double precision :: toxz
double precision :: toyy
double precision :: toyz
double precision :: tozz
double precision :: qcx
double precision :: qcy
double precision :: qcz
double precision :: tm1
double precision :: tm2
double precision :: tm3
double precision :: tm4
double precision :: tm5
double precision :: tm6
double precision :: tm7
double precision :: tm8
double precision :: tm9
double precision :: tm10
double precision :: tm11
double precision :: tm12
double precision :: tm13
integer :: ncin
integer :: mnc
integer :: ncbd
integer :: mnr
double precision :: xnr
double precision :: ynr
double precision :: znr
double precision :: bceqt
double precision :: rpi
double precision :: rti
double precision :: d0x
double precision :: d0y
double precision :: d0z
double precision :: qtx
double precision :: qty
double precision :: qtz
double precision :: pres
double precision :: tp
double precision :: rod
double precision :: roud
double precision :: rovd
double precision :: rowd
double precision :: roed
double precision :: pression
double precision :: ztemp
double precision :: cson
double precision :: cvi
double precision :: cvj
double precision :: cvk
double precision :: cmui1
double precision :: cmui2
double precision :: cmuj1
double precision :: cmuj2
double precision :: cmuk1
double precision :: cmuk2
integer :: icyexpl
integer :: img
integer :: iter
integer :: itypdf
integer :: ityprk
integer :: l
integer :: lm
integer :: mcyresi
integer :: mcysave
integer :: mfc
integer :: mfr
integer :: mglp
integer :: nc
integer :: ndeb
integer :: nfin
integer :: ngx
!
!-----------------------------------------------------------------------
!
      integer idcyc,ndcyc,icyc,ncyc,icycle,ncycl
      double precision nxn,nyn,nzn,mu,mut
!
      dimension mnpar(ip12),fgam(ip42),utau(ip42)
      dimension x(ip21),y(ip21),z(ip21)
      dimension dt(ip11),vol(ip11),r(ip11),pression(ip11),ztemp(ip11),cson(ip11)
      dimension v(ip11,ip60),d_volt(ip11,ip60),u_duv(ip11,ip60), &
             ff_du(ip11,ip60),u0(ip11,ip60)
      dimension vdual(ip11,ip60),vdual1(ip11,ip60),vdual2(ip11,ip60),ptdual(ip11,ip60)
      dimension sn(ip31*ndir)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
      dimension bceqt(ip41,neqt)
      dimension rpi(ip40),rti(ip40),pres(ip40),tp(ip40), &
                d0x(ip40),d0y(ip40),d0z(ip40), &
                qtx(ip40),qty(ip40),qtz(ip40),roed(ip40), &
                rod(ip40),roud(ip40),rovd(ip40),rowd(ip40)
      dimension xnr(ip44),ynr(ip44),znr(ip44),mnr(ip44)
      dimension mu(ip12),mut(ip12),toxx(ip12),toxy(ip12),toxz(ip12), &
                toyy(ip12),toyz(ip12),tozz(ip12),qcx(ip12),qcy(ip12), &
                qcz(ip12),dist(ip12)
      dimension cfke(ip13),ncin(ip41),mnc(ip43)
      dimension tn1(ip00),tn2(ip00),tn3(ip00),tn4(ip00),tn5(ip00), &
                tn6(ip00),tn7(ip00),tn8(ip00),tn9(ip00),tn10(ip00)
      dimension tm1(ip40),tm2(ip40),tm3(ip40),tm4(ip40),tm5(ip40), &
                tm6(ip40),tm7(ip40),tm8(ip40),tm9(ip40),tm10(ip40), &
                tm11(ip40),tm12(ip40),tm13(ip40)
      dimension cvi(ip21),cvj(ip21),cvk(ip21), &
                cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
!     niveau de grille courant --> mg
      mglp = mg+1
!
!-----------------------------------------------------------------------
!      Initialisation dual
!-----------------------------------------------------------------------
        if(kfmg.eq.3) then
          call sch_duin(v,1,ptdual,vdual,vdual1,vdual2)
        else
          niter = 1
          ndcyc = 0
!         test stockage pour calculs cavitants
          do l=1,lzx
           ndeb=npc(l)+1
           nfin=npc(l)+nnc(l)
           do nc=ndeb,nfin
            vdual2(nc,1)=0.
            vdual2(nc,2)=0.
            vdual2(nc,3)=0.
            vdual2(nc,4)=0.
            vdual2(nc,5)=0.
           enddo
         enddo
        endif
!
!*************************************************************************
!                                                         
!     INTEGRATION TEMPORELLE - boucle en temps physique                  
!                                                       
!*************************************************************************
!
      do icycle=1,ncycl
        if(kfmg.eq.3) then
          icyc=0
          ncyc=0
          ngx=lgx
        elseif(kfmg.eq.2) then
          icyc=idcyc
          ncyc=ndcyc
          ngx=mg
        elseif(kfmg.lt.2) then
          icyc=idcyc
          ncyc=ndcyc
          ngx=lgx
        endif
        idcyc = idcyc + 1
        ndcyc = icycle
        kdualto = 0
!
!-----------------------------------------------------------------------
!     boucle en temps dual
!-----------------------------------------------------------------------
!
       do iter=1,niter
        icyc=icyc+1
        ncyc=icycle+iter-1
        mcyresi=mod(ndcyc,ncyresi)
        mcysave=mod(ndcyc,ncysave)
        icyexpl=mod(idcyc,ncyexpl)
!
!-----------------------------------------------------------------
!     boucle sur niveau de grille
!-----------------------------------------------------------------
!
       do img = mg,ngx
!
        if(img.gt.mg) then     ! Grille grossiere en MG
         itypdf = 0
         if(img.gt.mglp) itypdf = 1
            call smg_flu( &
                 img-1,itypdf, &
                 ff_du, &
                 icyc,ncyc,icycle, &
                 x,y,z,r,nxn,nyn,nzn, &
                 sn, &
                 vol, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                 u_duv,mu,mut,dist, &
                 v,d_volt, &
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
!-------------------------------------------------------------------------
!     transfert des residus grille fine -> grille grossiere
!-------------------------------------------------------------------------
            call smg_fcr( &
                 img-1,img, &
                 ff_du)
!
!-------------------------------------------------------------------------
!     transfert des variables grille fine -> grille grossiere
!-------------------------------------------------------------------------
            call smg_fcv( &
                 img-1,img, &
                 vol,v)
      endif                            ! Fin grille grossiere MG
!
!-------------------------------------------------------------------------
!     stockage valeur initiale
!-------------------------------------------------------------------------
      do l = 1,lzx
       lm = l+(img-1)*lz
       ndeb=npc(lm)+1
       nfin=npc(lm)+nnc(lm)
!
       do nc=ndeb,nfin
        u0(nc,1)=v(nc,1)
        u0(nc,2)=v(nc,2)
        u0(nc,3)=v(nc,3)
        u0(nc,4)=v(nc,4)
        u0(nc,5)=v(nc,5)
       enddo
!
       if (equat(6:7).eq.'ke') then
        do nc=ndeb,nfin
         u0(nc,6)=v(nc,6)
         u0(nc,7)=v(nc,7)
        enddo
       endif
      enddo
!
      ityprk = 1
      if(img.eq.mg) ityprk = 0
!
!-------------------------------------------------------------------------
!       application du schema - avance d'un pas de temps
!-------------------------------------------------------------------------
!
        call atsch_num( &
                 img,ityprk, &
                 icyc,ncyc,idcyc,icycle, &
                 mg, &
                 x,y,z,r,nxn,nyn,nzn, &
                 sn, &
                 vol, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                 u_duv,mu,mut,dist,cfke, &
                 mnpar,fgam,utau, &
                 v,dt,d_volt, &
                 ptdual, &
                 ff_du, &
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
      enddo                       ! Fin boucle sur img (descente)
!
      img=mg
!
      if((kfmg.eq.3).and.(kdualto.eq.0)) then
        call dua_resro( &
                 icyc,ncyc,img, &
                 u0,v,dt)
!
        if(resite.le.tol) kdualto = 1
      endif
!
!-------------------------------------------------------------------------
!  calcul des residus et eventuellement calculs utilisateur
!-------------------------------------------------------------------------
      if(mcyresi.eq.0) then
        if(img.eq.mg) then
          call rscpsv( &
                 img, &
                 u0,v,dt, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8, &
                 icyc,ncyc,ncycl, &
                 x,y,z,utau)
        endif
      endif
!
      if(kfmg.ne.2) then
       if(mglp.le.lgx) then
        do img = lgx,mglp,-1
!-------------------------------------------------------------------------
!     corrections (level=img)
!-------------------------------------------------------------------------
          do l = 1,lzx
           lm = l+(img-1)*lz
           ndeb=npc(lm)+1
           nfin=npc(lm)+nnc(lm)
           do nc=ndeb,nfin
            ff_du(nc,1)=v(nc,1)-u0(nc,1)
            ff_du(nc,2)=v(nc,2)-u0(nc,2)
            ff_du(nc,3)=v(nc,3)-u0(nc,3)
            ff_du(nc,4)=v(nc,4)-u0(nc,4)
            ff_du(nc,5)=v(nc,5)-u0(nc,5)
           enddo
          enddo
!
          if(mtcx.gt.0) then
           do mfc=1,mtcx
            lbd(mfc)=nfbc(mfc)+(img-1)*mtb
           enddo
           nbd=mtcx
!           call rfvc(ff_du,ncbd,mnc)
           call rfvc( & 
                ff_du,ncbd,mnc, &
                pression,ztemp,cson)
          endif
!
          if (mtrx.gt.0) then
           do mfr=1,mtrx
            lbd(mfr)=nfbr(mfr)+(img-1)*mtb
           enddo
           nbd=mtrx
           call rfvr( &
                 ff_du, &
                 ncbd,mnr,xnr,ynr,znr)
          endif
!-------------------------------------------------------------------------
!     cell-to-node projection (level=img)
!-------------------------------------------------------------------------
            call smg_cn( &
                 img, &
                 vol,d_volt(1,1), &
                 ff_du,u_duv)
!-------------------------------------------------------------------------
!      Node-to-cell interpolation (level=img-1)
!-------------------------------------------------------------------------
            call smg_cf( &
                 img,img-1, &
                 vol, &
                 u_duv,ff_du)
!-------------------------------------------------------------------------
!      Mise a jour des corrections (level=img-1)
!-------------------------------------------------------------------------
            call smg_upc( &
                 img-1, &
                 ff_du,v)
           enddo                       ! fin boucle img
          endif
         endif                         ! fin test sur kfmg
!
!      kdualns=0 : DTS + Euler
!      kdualns=1 : DTS + NS turbulent (eq. transport) ordre 2 en temps
!      kdualns=2 : DTS + NS turbulent (eq. transport) ordre 3 en temps
!
        if((kfmg.eq.3).and.(kdualto.eq.1).and.(kdualns.eq.0)) exit
        if((kfmg.eq.3).and.(kdualto.eq.1).and.(kdualns.eq.1)) exit
        if((kfmg.eq.3).and.(kdualto.eq.1).and.(kdualns.eq.2)) exit
!
       enddo   ! Fin boucle temps dual
!
!-------------------------------------------------------------------------
!      Sorties et postraitement
!-------------------------------------------------------------------------
!
      if(mg.eq.1) then
       if(icyexpl.eq.0) then
!
!     Fabrication d'un champs verifiant les conditions aux limites
!     pour exploitation eventuelle par l'utilisateur
!
            call rbc( &
                 ncin,nxn,nyn,nzn,ncbd, &
                 sn,vol,v,mut, &
                 bceqt, &
                 rpi,rti,d0x,d0y,d0z,qtx,qty,qtz,x,y,z,omg, &
                 pres,tp,rod,roud,rovd,rowd,roed, &
                 icyc, &
                 mnr,xnr,ynr,znr,mnc, &
                 tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10,tm11, &
                 tm12,tm13,pression,ztemp,cson)
        endif
!
         if(equat(1:2).eq.'ns') then ! calcul Navier-Stokes

            call utitfr(  &
                x,y,z,v,ncbd,nxn,nyn,nzn,idcyc,  &
                toxx,toxy,toxz,toyy,toyz,tozz,  &
                pression)

           if((kfmg.eq.3).and.(lsortie.eq.1)) then
             do l=1,lzx  !sigma / moyenne temporelle / pression
              call sortietest(   &
                 icycle,ncycl,idcyc, &
                 vdual2,dist,vol,mut,mu, &
                 x,y,z,l,v,pression,ztemp)
             enddo
!
!          sortie instationnaire densite-pression
            if(mod(idcyc,nfreq).eq.0) then
             do l=1,lzx
              call sortieplot(x,y,z,l,v,pression,cson)
!              call sortieplot(x,y,z,l,v,pression,dist,mu,mut)
             enddo
            endif
           endif
!
          else ! calcul Euler
!
!            call utit(x,y,z,v,ncbd,nxn,nyn,nzn,idcyc)

           if((kfmg.eq.3).and.(lsortie.eq.1)) then
!          sortie instationnaire densite-pression
            if(mod(idcyc,nfreq).eq.0) then
             do l=1,lzx
              call sortieplot(x,y,z,l,v,pression,cson)
              enddo
             endif
           endif

        endif
      endif     ! Fin test sur mg
!
!-------------------------------------------------------------------------
!        sauvegarde du fichier resultats
!-------------------------------------------------------------------------
      if(((mcysave.eq.0).or.(ndcyc.eq.ncycl)).and.(img.eq.1)) then
!      if((ndcyc.eq.ncycl).and.(img.eq.1)) then
!
       do l=1,lzx
            call svfw( &
                 l,v,mut,utau, &
                 'cccc',ncin,ncbd, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8)
!
          if(kfmg.lt.3) then  
!         sorties tecplot
            call sortieplot2(    &
               x,y,z,l,v,dist, &
               mu,mut,toxy, &
               pression,cson,ztemp)
!            call sortieplot3(    &
!               x,y,z,l,v,dist, &
!               mu,mut,toxy, &
!               pression,cson,ztemp)
          endif
       enddo
      endif
!-------------------------------------------------------------------------
!      Mise a jour tableaux vdual, vdual1, vdual2, ptdual
!-------------------------------------------------------------------------
       if(kfmg.eq.3) then
         if((kdualns.eq.0).or.(kdualns.eq.1)) then
!         if(kdualns.eq.1) then
           call sch_duup( &
                  sn,vol,tn3, &
                  v,ptdual,vdual,vdual1, &
                  mu,mut,dist,tn4, &
                  ncyc,1, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
           icyc = 0
         elseif(kdualns.eq.2) then
!         elseif((kdualns.eq.0).or.(kdualns.eq.2)) then
            call sch_duup2( &
                  sn,vol,tn3, &
                  v,ptdual,vdual,vdual1,vdual2, &
                  mu,mut,dist,tn4, &
                  ncyc,1, &
                  cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
           icyc = 0
         endif
       endif
!
      enddo    ! fin boucle icycle sur temps physique
!
      return
      end subroutine
end module
