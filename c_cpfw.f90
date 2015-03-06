      subroutine c_cpfw( &
                 mot,imot,nmot, &
                 ncyc, &
                 x,y,z,r,exs1,exs2,nxn,nyn,nzn, &
                 sn, &
                 vol, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                 mu,mut,dist,cfke, &
                 mnpar,fgam,utau, &
                 v,dt, &
                 ptdual,vdual,vdual1,vdual2, &
                 tnte1,tnte2,tnte3,tnte4, &
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
!     ACT
!_A    Realisation de l'action cpfw.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use sortiefichier
!
!-----------------------------------------------------------------------
!
      character *32 mot(nmx)
      real nxn,nyn,nzn,mu,mut
!
      dimension imot(nmx)
      dimension x(ip21),y(ip21),z(ip21)
      dimension v(ip11,ip60)
      dimension dt(ip11),vol(ip11),r(ip11),pression(ip11),ztemp(ip11),cson(ip11)
      dimension sn(ip31*ndir)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
      dimension bceqt(ip41,neqt)
      dimension rpi(ip40),rti(ip40),d0x(ip40),d0y(ip40),d0z(ip40)
      dimension qtx(ip40),qty(ip40),qtz(ip40),pres(ip40),tp(ip40)
      dimension rod(ip40),roud(ip40),rovd(ip40),rowd(ip40),roed(ip40)
      dimension mu(ip12),mut(ip12),toxx(ip12),toxy(ip12),toxz(ip12), &
                toyy(ip12),toyz(ip12),tozz(ip12),qcx(ip12),qcy(ip12), &
                qcz(ip12),dist(ip12),mnpar(ip12)
      dimension fgam(ip42),utau(ip42),cfke(ip13),mnc(ip43),ncin(ip41)
      dimension tnte1(ip11,ip60),tnte2(ip11,ip60), &
                tnte3(ip11,ip60),tnte4(ip11,ip60)
      dimension tn1(ip00),tn2(ip00),tn3(ip00),tn4(ip00),tn5(ip00), &
                tn6(ip00),tn7(ip00),tn8 (ip00),tn9(ip00),tn10(ip00)
      dimension tm1(ip40),tm2(ip40),tm3(ip40),tm4(ip40),tm5(ip40), &
                tm6(ip40),tm7(ip40),tm8(ip40),tm9(ip40),tm10(ip40), &
                tm11(ip40),tm12(ip40),tm13(ip40)
      dimension xnr(ip44),ynr(ip44),znr(ip44),mnr(ip44)
      dimension vdual(ip11,ip60),vdual1(ip11,ip60),vdual2(ip11,ip60),ptdual(ip11,ip60)
      dimension cvi(ip21),cvj(ip21),cvk(ip21), &
                cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
                cmuk1(ip21),cmuk2(ip21)
!
      call tcmd_cpfw(mot,imot,nmot)
!
      if(kimp.ge.1) then
            call b1_cpfw
      endif
!
      call cpfw( &
                 ncyc, &
                 x,y,z,r,exs1,exs2,nxn,nyn,nzn, &
                 sn, &
                 vol, &
                 tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
                 mu,mut,dist,cfke, &
                 mnpar,fgam,utau, &
                 v,dt, &
                 ptdual,vdual,vdual1,vdual2, &
                 tnte1,tnte2,tnte3,tnte4, &
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
      return
      end
