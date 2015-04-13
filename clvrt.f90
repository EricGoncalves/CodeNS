module mod_clvrt
  implicit none
contains
  subroutine clvrt( &
       mfb, &
       rod,roud,rovd,rowd,roed, &
       nxn,nyn,nzn,ncbd,v, &
       mmb,mpb,mpn, &
       l,x,y,z, &
       pression,temp,cson)
!
!***********************************************************************
!
!    ACT
!     condition limite utilisable en transsonique pour restituer les
!     conditions a l'infini sur un maillage limite.
!
!***********************************************************************
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    use schemanum
    implicit none
    integer          ::         id,        jd,        kd,         l,         m
    integer          ::         mb,       mfb,  mmb(mtt),        mn,  mpb(mtt)
    integer          ::   mpn(mtt),        mt,       n0c,       n0n,ncbd(ip41)
    integer          ::        nci,      ncij,     ncijk,      ncik,       ncj
    integer          ::       ncjk,       nck,       nid,      nijd,       njd
    integer          ::         nl,        nn
    double precision ::             am,            ap,            b0,            bs,    cson(ip11)
    double precision ::            dvx,           dvz,            dx,            dz,          eps0
    double precision ::           epsm,          epsp,          gami,         gm1sg,         gsgm1
    double precision ::      nxn(ip42),     nyn(ip42),     nzn(ip42),             p,           pd0
    double precision ::            pdd,            pi,pression(ip11),            ps,         qcs4p
    double precision ::            qd0,           qd2,          qd20,            qn,           qnd
    double precision ::            qns,           qtx,          qtxd,          qtxs,           qty
    double precision ::           qtyd,          qtys,           qtz,          qtzd,          qtzs
    double precision ::            qxd,          qxd0,           qxs,           qyd,          qyd0
    double precision ::            qys,           qzd,          qzd0,           qzs,            ro
    double precision ::            ro0,          roc0,     rod(ip40),          rod0,          rodd
    double precision ::     roed(ip40),         roqn0,           ros,    roud(ip40),    rovd(ip40)
    double precision ::     rowd(ip40),    temp(ip11),  v(ip11,ip60),        vrtalr,        vrtang
    double precision ::         vrtbet,        vrtdel,        vrtden,        vrtdis,        vrtfac
    double precision ::         vrtsin,          vrtx,          vrtz,       x(ip21),       y(ip21)
    double precision ::        z(ip21)
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!

!
    n0n=npn(l)
    n0c=npc(l)
!
    nid = id2(l)-id1(l)+1
    njd = jd2(l)-jd1(l)+1
    nijd = nid*njd
!
    nci = inc(1,0,0)
    ncj = inc(0,1,0)
    nck = inc(0,0,1)
    ncij = inc(1,1,0)
    ncik = inc(1,0,1)
    ncjk = inc(0,1,1)
    ncijk= inc(1,1,1)
!
    mt=mmb(mfb)
!
    gm1sg=(gam-1.)/gam
    gami=1./gam
    gsgm1=gam/(gam-1.)
!
    pi=4.*atan(1.)
    vrtalr=pi*vrtalp/180.
    if(vrtmac.ge.1.)then
       print*,'!!!clvrt: condition vort interdite en supersonique'
       stop
    else
       vrtbet=sqrt(1.-vrtmac*vrtmac)
    endif
!
!$OMP SIMD
    do m=1,mt
!
       qxd0=roud(m)/rod(m)
       qyd0=rovd(m)/rod(m)
       qzd0=rowd(m)/rod(m)
       qd20=qxd0**2+qyd0**2+qzd0**2
       pd0=gam1*(roed(m)-.5*rod(m)*qd20)
       qd0=sqrt(qd20)
       rod0=rod(m)
       qcs4p=qd0*vrtlre*.25/pi
!
       mb=mpb(mfb)+m
       mn=mpn(mfb)+m
       nl=ncbd(mb)
       nn=nl-n0c+n0n
!
       vrtx=0.125*( x (nn     )+x (nn+nci  ) &
            +x (nn+ncj )+x (nn+ncij ) &
            +x (nn+nck )+x (nn+ncik ) &
            +x (nn+ncjk)+x (nn+ncijk) )
       vrtz=0.125*( z (nn     )+z (nn+nci  ) &
            +z (nn+ncj )+z (nn+ncij ) &
            +z (nn+nck )+z (nn+ncik ) &
            +z (nn+ncjk)+z (nn+ncijk) )
!
       dx=vrtx-vrtxre
       dz=vrtz-vrtzre
       vrtdis=sqrt(dx*dx+dz*dz)
       vrtdel=atan2(dz,dx)-vrtalr
       vrtsin=sin(vrtdel)
       vrtden=1.-vrtmac*vrtmac*vrtsin*vrtsin
       vrtfac=qcs4p*vrtbet/vrtden/vrtdis
       vrtang=atan2(dz,dx)
       dvx=vrtcz*vrtfac*sin(vrtang)
       dvz=-vrtcz*vrtfac*cos(vrtang)
!
       ro0 =v(nl,1)
       roc0=gam5*(ro0*v(nl,5) &
            -.5*(v(nl,2)**2+v(nl,3)**2+v(nl,4)**2))
       roc0=sqrt(roc0)
       roqn0=v(nl,2)*nxn(mn)+v(nl,3)*nyn(mn)+v(nl,4)*nzn(mn)
       epsm=.5+sign(.5, roc0-roqn0)
       eps0=.5+sign(.5,-roqn0)
       epsp=.5+sign(.5,-roc0-roqn0)
!
       qxd=roud(m)/rod(m)+dvx
       qyd=rovd(m)/rod(m)
       qzd=rowd(m)/rod(m)+dvz
!
       qd2=qxd*qxd+qyd*qyd+qzd*qzd
       pdd=(pd0**gm1sg+gm1sg*rod0*(qd20-qd2)/(2.*pd0**gami))**gsgm1
       rodd=rod0*(pdd/pd0)**gami
!
       qnd =qxd*nxn(mn)+qyd*nyn(mn)+qzd*nzn(mn)
       qtxd=qxd-qnd*nxn(mn)
       qtyd=qyd-qnd*nyn(mn)
       qtzd=qzd-qnd*nzn(mn)
!
       ros =v(nl,1)
       qxs =v(nl,2)/ros
       qys =v(nl,3)/ros
       qzs =v(nl,4)/ros
       ps  =gam1*(v(nl,5)-.5*ros*(qxs**2+qys**2+qzs**2))
       qns =qxs*nxn(mn)+qys*nyn(mn)+qzs*nzn(mn)
       qtxs=qxs-qns*nxn(mn)
       qtys=qys-qns*nyn(mn)
       qtzs=qzs-qns*nzn(mn)
!
       qtx=eps0*qtxs         +(1.-eps0)*qtxd
       qty=eps0*qtys         +(1.-eps0)*qtyd
       qtz=eps0*qtzs         +(1.-eps0)*qtzd
       am =epsm*(ps-roc0*qns)+(1.-epsm)*(pdd-roc0*qnd)
       ap =epsp*(ps+roc0*qns)+(1.-epsp)*(pdd+roc0*qnd)
       qn  =(ap-am)*.5/roc0
       p   =(ap+am)*.5
       bs  =(p-ps  )*ro0**2/roc0**2+ros
       b0  =(p-pdd)*ro0**2/roc0**2+rodd
       ro  =eps0*bs           +(1.-eps0)*b0
!
       v(nl,1)=ro
       v(nl,2)=ro*(qtx+qn*nxn(mn))
       v(nl,3)=ro*(qty+qn*nyn(mn))
       v(nl,4)=ro*(qtz+qn*nzn(mn))
       v(nl,5)=p*gam4+0.5*(v(nl,2)**2+v(nl,3)**2+v(nl,4)**2)/ro
       pression(nl)=p
       temp(nl)=gam*pression(nl)/ro
       cson(nl)=sqrt(temp(nl))
    enddo
!
!$OMP END MASTER
  contains
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine clvrt

end module mod_clvrt
