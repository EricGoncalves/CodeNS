module mod_c_cpfw
  implicit none
contains
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
    use mod_b1_cpfw

    use mod_tcmd_cpfw
    use mod_cpfw
    implicit none
    integer          ::   imot(nmx),  mnc(ip43),mnpar(ip12),  mnr(ip44), ncbd(ip41)
    integer          ::  ncin(ip41),       ncyc,       nmot
    double precision ::  bceqt(ip41,neqt),       cfke(ip13),      cmui1(ip21),      cmui2(ip21),      cmuj1(ip21)
    double precision ::       cmuj2(ip21),      cmuk1(ip21),      cmuk2(ip21),       cson(ip11),        cvi(ip21)
    double precision ::         cvj(ip21),        cvk(ip21),        d0x(ip40),        d0y(ip40),        d0z(ip40)
    double precision ::        dist(ip12),         dt(ip11),             exs1,             exs2,       fgam(ip42)
    double precision ::          mu(ip12),        mut(ip12),        nxn(ip42),        nyn(ip42),        nzn(ip42)
    double precision ::        pres(ip40),   pression(ip11),ptdual(ip11,ip60),        qcx(ip12),        qcy(ip12)
    double precision ::         qcz(ip12),        qtx(ip40),        qty(ip40),        qtz(ip40),          r(ip11)
    double precision ::         rod(ip40),       roed(ip40),       roud(ip40),       rovd(ip40),       rowd(ip40)
    double precision ::         rpi(ip40),        rti(ip40),    sn(ip31*ndir),        tm1(ip40),       tm10(ip40)
    double precision ::        tm11(ip40),       tm12(ip40),       tm13(ip40),        tm2(ip40),        tm3(ip40)
    double precision ::         tm4(ip40),        tm5(ip40),        tm6(ip40),        tm7(ip40),        tm8(ip40)
    double precision ::         tm9(ip40),        tn1(ip00),       tn10(ip00),        tn2(ip00),        tn3(ip00)
    double precision ::         tn4(ip00),        tn5(ip00),        tn6(ip00),        tn7(ip00),        tn8(ip00)
    double precision ::         tn9(ip00), tnte1(ip11,ip60), tnte2(ip11,ip60), tnte3(ip11,ip60), tnte4(ip11,ip60)
    double precision ::        toxx(ip12),       toxy(ip12),       toxz(ip12),       toyy(ip12),       toyz(ip12)
    double precision ::        tozz(ip12),         tp(ip40),       utau(ip42),     v(ip11,ip60), vdual(ip11,ip60)
    double precision :: vdual1(ip11,ip60),vdual2(ip11,ip60),        vol(ip11),          x(ip21),        xnr(ip44)
    double precision ::           y(ip21),        ynr(ip44),          z(ip21),        znr(ip44),      ztemp(ip11)
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
!
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
  end subroutine c_cpfw
end module mod_c_cpfw
