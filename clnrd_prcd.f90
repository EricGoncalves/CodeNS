module mod_clnrd_prcd
  implicit none
contains
  subroutine clnrd_prcd( &
       mfb, &
       rod,roud,rovd,rowd,roed, &
       nxn,nyn,nzn,ncbd,v, &
       mmb,mpb,mpn,l, &
       pression,temp,cson)
!
!***********************************************************************
!
!_DA           aout 2002 - Eric GONCALVES / SINUMEF
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement d'une condition de non reflexion.
!_A    Utilisation des cinq relations de compatibilite discretisees
!_A    combinees avec la donnee d'un etat de reference (utnrd.f).
!_A      Dut=0
!_A      Dut=0
!_A      Dp - c2*Drho = 0
!_A     (lambda- -u)*Dp + rho*beta2*c2*Du = 0
!_A     (lambda+ -u)*Dp + rho*beta2*c2*Du = 0
!_A    Normales interieures.
!_A    Preconditionnement basse vitesse de Turkel.
!
!     INP
!_I    mfb        : arg int              ; numero de frontiere
!_I    rod        : arg real(ip40      ) ; masse volumique donnee
!_I    roud       : arg real(ip40      ) ; ro*u donne, ro masse volumique
!_I                                        et u composante en x de la vitesse
!_I    rovd       : arg real(ip40      ) ; ro*v donne, ro masse volumique
!_I                                        et v composante en y de la vitesse
!_I    rowd       : arg real(ip40      ) ; ro*w donne, ro masse volumique
!_I                                        et w composante en z de la vitesse
!_I    roed       : arg real(ip40      ) ; ro*e donne, ro masse volumique
!_I                                        et e energie interne
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!
!     COM
!_C    Notation s      : valeurs issues du schema.
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use proprieteflu
    use schemanum
    use definition
    implicit none
    integer          ::          l,         m,        mb,       mfb,  mmb(mtt)
    integer          ::         mn,  mpb(mtt),  mpn(mtt),        mt,       n0c
    integer          ::        n0n,ncbd(ip41),        nl
    double precision ::             a2,           alm,           alp,            am,            ap
    double precision ::             b0,         beta2,            bs,    cson(ip11),          eps0
    double precision ::           epsm,          epsp,     nxn(ip42),     nyn(ip42),     nzn(ip42)
    double precision ::              p,            pd,pression(ip11),            ps,            q2
    double precision ::           qinf,            qn,           qnd,           qns,           qtx
    double precision ::           qtxd,          qtxs,           qty,          qtyd,          qtys
    double precision ::            qtz,          qtzd,          qtzs,           qxd,           qxs
    double precision ::            qyd,           qys,           qzd,           qzs,           rho
    double precision ::           rhos,           rmd,           rms,          rocs,     rod(ip40)
    double precision ::     roed(ip40),    roud(ip40),    rovd(ip40),    rowd(ip40),           rpd
    double precision ::            rps,    temp(ip11),  v(ip11,ip60)
!
!-----------------------------------------------------------------------
!
!
!
    n0n=npn(l)
    n0c=npc(l)
    mt=mmb(mfb)
    qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
!
    do m=1,mt
       mb=mpb(mfb)+m
       mn=mpn(mfb)+m
       nl=ncbd(mb)
!      etat schema
       rhos=v(nl,1)
       qxs=v(nl,2)/rhos
       qys=v(nl,3)/rhos
       qzs=v(nl,4)/rhos
       ps=gam1*(v(nl,5)-0.5*rhos*(qxs**2+qys**2+qzs**2)-pinfl)
       qns=qxs*nxn(mn)+qys*nyn(mn)+qzs*nzn(mn)
       rocs=sqrt(gam*ps*rhos)
       q2=qxs**2+qys**2+qzs**2
       a2=gam*ps/rhos
       beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
       alp=0.5*((1.+beta2)*abs(qns) + sqrt(((1.-beta2)*qns)**2+4.*beta2*a2))
       alm=0.5*((1.+beta2)*abs(qns) - sqrt(((1.-beta2)*qns)**2+4.*beta2*a2))
       eps0=.5+sign(.5,-qns)
       epsp=.5+sign(.5,-alp)
       epsm=.5+sign(.5,-alm)
       qtxs=qxs-qns*nxn(mn)
       qtys=qys-qns*nyn(mn)
       qtzs=qzs-qns*nzn(mn)
!      etat de reference (valeur prise a t=0)
       qxd=roud(m)/rod(m)
       qyd=rovd(m)/rod(m)
       qzd=rowd(m)/rod(m)
       pd=gam1*(roed(m)-0.5*rod(m)*(qxd**2+qyd**2+qzd**2)-pinfl)
       qnd=qxd*nxn(mn)+qyd*nyn(mn)+qzd*nzn(mn)
       qtxd=qxd-qnd*nxn(mn)
       qtyd=qyd-qnd*nyn(mn)
       qtzd=qzd-qnd*nzn(mn)
!
       qtx=eps0*qtxs   +(1.-eps0)*qtxd
       qty=eps0*qtys   +(1.-eps0)*qtyd
       qtz=eps0*qtzs   +(1.-eps0)*qtzd
       rms=(alm-qns)*ps+rhos*beta2*a2*qns
       rmd=(alm-qns)*pd+rhos*beta2*a2*qnd
       rps=(alp-qns)*ps+rhos*beta2*a2*qns
       rpd=(alp-qns)*pd+rhos*beta2*a2*qnd
       am=epsm*rms     +(1.-epsm)*rmd
       ap=epsp*rps     +(1.-epsp)*rpd
       qn=((alm-qns)*ap-(alp-qns)*am)/((alm-alp)*rhos*beta2*a2)
       p=(am-ap)/(alm-alp)
       bs=(p-ps)*rhos**2/rocs**2+rhos
       b0=(p-pd)*rhos**2/rocs**2+rod(m)
       rho=eps0*bs     +(1.-eps0)*b0
!
       v(nl,1)=rho
       v(nl,2)=rho*(qtx+qn*nxn(mn))
       v(nl,3)=rho*(qty+qn*nyn(mn))
       v(nl,4)=rho*(qtz+qn*nzn(mn))
       v(nl,5)=p/gam1+pinfl+0.5*(v(nl,2)**2+v(nl,3)**2+v(nl,4)**2)/rho
!
       pression(nl)=p
       temp(nl)=gam*pression(nl)/rho
       cson(nl)=sqrt(temp(nl))
    enddo
!
    return
  end subroutine clnrd_prcd
end module mod_clnrd_prcd
