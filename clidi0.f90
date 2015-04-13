module mod_clidi0
  implicit none
contains
  subroutine clidi0( &
       mfb,l,rpi,rti, &
       nxn,nyn,nzn,ncbd,v, &
       y,z, &
       mmb,mpb,mpn, &
       usdn2,roc0,am0,qn,p,resi,ro,rovn0,un,vy0,vz0, &
       ym,zm,pression,temp,cson)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement des conditions aux limites :
!_A     (direction de la vitesse absolue)/dt=0,
!_A     pression d'arret imposee ,  pi=pa1*rpi
!_A     enthalpie d'arret imposee , hi=ha1*rti
!_A    Utilisation d une relation de compatibilite en ro,p,qx,qy,qz.
!_A    Normales interieures.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    use definition
    implicit none
    integer          ::         id,        jd,        kd,         l,         m
    integer          ::         mb,       mfb,  mmb(mtt),        mn,  mpb(mtt)
    integer          ::   mpn(mtt),        mt,       n0c,       n0n,        nc
    integer          :: ncbd(ip41),       nci,      ncij,     ncijk,      ncik
    integer          ::        ncj,      ncjk,       nck,       nid,      nijd
    integer          ::       nitn,       njd,        nn
    double precision ::      am0(ip40),             b,          coef,    cson(ip11),            df
    double precision ::            dqn,           eps,             f,         gam2t,          gam6
    double precision ::           gamt,     nxn(ip42),     nyn(ip42),     nzn(ip42),       p(ip40)
    double precision :: pression(ip11),            ps,      qn(ip40),           qns,           qxs
    double precision ::            qys,           qzs,    resi(ip40),        residu,      ro(ip40)
    double precision ::     roc0(ip40),         roq20,         roqn0,   rovn0(ip40),         rovv0
    double precision ::      rpi(ip40),     rti(ip40),    temp(ip11),      un(ip40),   usdn2(ip40)
    double precision ::   v(ip11,ip60),            vn,     vy0(ip40),     vz0(ip40),       y(ip21)
    double precision ::       ym(ip40),       z(ip21),      zm(ip40)
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!

    eps=0.0000001
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
    gam6=gam/gam1
    gamt=gam/aa1**2
    gam2t=gam2/aa1**2
!
!$OMP SIMD
    do m=1,mt
       mb=mpb(mfb)+m
       mn=mpn(mfb)+m
       nc=ncbd(mb)
       nn=nc-n0c+n0n
!
       ym(m) = 0.125*( y (nn     )+y (nn+nci  ) &
            +y (nn+ncj )+y (nn+ncij ) &
            +y (nn+nck )+y (nn+ncik ) &
            +y (nn+ncjk)+y (nn+ncijk) )
       zm(m) = 0.125*( z (nn     )+z (nn+nci  ) &
            +z (nn+ncj )+z (nn+ncij ) &
            +z (nn+nck )+z (nn+ncik ) &
            +z (nn+ncjk)+z (nn+ncijk) )
!
       roqn0=v(nc,2)*nxn(mn)+v(nc,3)*nyn(mn)+v(nc,4)*nzn(mn)
       roq20=v(nc,2)**2+v(nc,3)**2+v(nc,4)**2
       qn(m)=roqn0/v(nc,1)
!
       vy0(m)=v(nc,3)-omg*zm(m)*v(nc,1)
       vz0(m)=v(nc,4)+omg*ym(m)*v(nc,1)
       rovn0(m)=v(nc,2)*nxn(mn)+vy0(m)*nyn(mn)+vz0(m)*nzn(mn)
       rovv0=v(nc,2)**2+vy0(m)**2+vz0(m)**2
       usdn2(m)=rovv0/rovn0(m)**2
       un(m)=omg*(ym(m)*nzn(mn)-zm(m)*nyn(mn))
!
       roc0(m)=gam5*(v(nc,1)*v(nc,5)-.5*roq20)
       roc0(m)=-sqrt(roc0(m))
!
       qxs=v(nc,2)/v(nc,1)
       qys=v(nc,3)/v(nc,1)
       qzs=v(nc,4)/v(nc,1)
       qns=qxs*nxn(mn)+qys*nyn(mn)+qzs*nzn(mn)
       ps=gam1*(v(nc,5)-0.5*(qxs*v(nc,2)+qys*v(nc,3)+qzs*v(nc,4))-pinfl)
       am0(m)=ps+roc0(m)*qns
    enddo
!
!   resolution de p+rcm0*qn-am0=0 par newton :
!
    nitn=0
    residu=1.
    do while( (residu.gt.eps).and.(nitn.lt.20) )
       residu=0.
!$OMP SIMD
       do m=1,mt
          mb=mpb(mfb)+m
!-------la notation b designe la quantite 1/(1+((gam-1)/2)m**2)
          b=1./(1.+gam2t/rti(m)*usdn2(m)*(qn(m)+un(m))**2)
          p(m)=pa1*rpi(m)*b**gam6
          ro(m)=gamt*p(m)/(b*rti(m))
          f=p(m)+roc0(m)*qn(m)-am0(m)
          df=roc0(m)-ro(m)*(qn(m)+un(m))*usdn2(m)
          dqn=-f/df
          resi(m)=abs(dqn)/qn(m)
          qn(m)=qn(m)+dqn

          residu=max(residu,resi(m))
       enddo
       nitn=nitn+1
    enddo
!
!$OMP SIMD
    do m=1,mt
       mb=mpb(mfb)+m
       nc=ncbd(mb)
!
       b=1./(1.+gam2t/rti(m)*usdn2(m)*(qn(m)+un(m))**2)
       p(m)=pa1*rpi(m)*b**gam6
       ro(m)=gamt*p(m)/(b*rti(m))
       vn=qn(m)+un(m)
       coef=vn/rovn0(m)
!
       v(nc,1)=ro(m)
       v(nc,2)=ro(m)*(coef*v(nc,2))
       v(nc,3)=ro(m)*(coef*vy0(m)+zm(m)*omg)
       v(nc,4)=ro(m)*(coef*vz0(m)-ym(m)*omg)
       v(nc,5)=p(m)/gam1+pinfl+0.5*(v(nc,2)**2+v(nc,3)**2+v(nc,4)**2)/ro(m)
!
       pression(nc)=p(m)
       temp(nc)=gam*p(m)/ro(m)
       cson(nc)=sqrt(temp(nc))
    enddo
!
!$OMP END MASTER
    return
  contains
    function    inc(id,jd,kd)
      implicit none
      integer          ::  id,inc, jd, kd
      inc=id+jd*nid+kd*nijd
    end function inc
  end subroutine clidi0
end module mod_clidi0
