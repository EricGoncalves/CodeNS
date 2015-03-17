module mod_clidd2_prcd
  implicit none
contains
  subroutine clidd2_prcd( &
       mfb,l,rpi,rti,d0x,d0y,d0z, &
       nxn,nyn,nzn,ncbd,v, &
       y,z, &
       mmb,mpb,mpn, &
       usdn2,roc0,am0,qn,p,resi,ro,un,usdn, &
       ym,zm,alm,pression,temp,cson)
!
!***********************************************************************
!
!_DA   avril 2003 - AUTEUR: Eric GONCALVES / LEGI
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement des conditions aux limites :
!_A     direction de la vitesse absolue imposee, parallele a (d0x,d0y,d0z),
!_A     pression d'arret imposee , pi=pa1*rpi
!_A     enthalpie d'arret imposee ,  hi=ha1*rti
!_A    Utilisation de la relation de compatibilite discretisee associee
!_A    a la caracteristique sortante (lambda-)
!_A     (lambda- -u)*Dp + rho*beta2*c2*Du = 0
!_A    Normales interieures.
!_A    Preconditionnement basse vitesse de Turkel.
!
!     INP
!_I    ip11       : arg int              ; dim, nbr max de cellules de tous les
!_I                                        dom (pts fictifs inclus)
!_I    ip21       : arg int              ; dim, nbr max de noeuds de tous les
!_I                                        dom (pts fictifs inclus)
!_I    ip40       : arg int              ; dim, nbr max de pts d'une front
!_I    ip41       : arg int              ; dim, nbr max de pts de ttes les front
!_I    ip42       : arg int              ; dim, nbr max de pts de ttes les front
!_I                                        a normales stockees
!_I    ip60       : arg int              ; dim, nbr max d'equations
!_I    mfb        : arg int              ; numero de frontiere
!_I    l          : arg int              ; numero de domaine
!_I    rpi        : arg real(ip40      ) ; pres d'arret/pres d'arret etat de
!_I                                        ref utilisateur a imposer
!_I    rti        : arg real(ip40      ) ; temp d'arret/temp d'arret etat de
!_I                                        ref utilisateur a imposer
!_I    d0x        : arg real(ip40      ) ; composante en x de la
!_I                                        direction de l'ecoulement a imposer
!_I    d0y        : arg real(ip40      ) ; composante en y de la
!_I                                        direction de l'ecoulement a imposer
!_I    d0z        : arg real(ip40      ) ; composante en z de la
!_I                                        direction de l'ecoulement a imposer
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    omg        : arg real             ; vitesse rotation du repere relatif
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    aa1        : com real             ; etat de reference utilisateur
!_I                                        adimensionne, vitesse du son d'arret
!_I    pa1        : com real             ; pression d'arret de l'etat
!_I                                        de reference utilisateur adimensionne
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    npn        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab tous noeuds
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!
!     I/O
!_/    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!
!     LOC
!_L    usdn2      : arg real(ip40      ) ; 1/(d.n)**2, d direction de la vitesse
!_L                                        et n normale a la front
!_L    roc0       : arg real(ip40      ) ; ro*c0, ro masse volumique,
!_L                                        c0 vitesse du son
!_L    am0        : arg real(ip40      ) ; ps+ro*c0*qn, ps pression statique
!_L                                        schema, ro masse volumique, c0 vit.
!_L                                        du son, qn vitesse normale schema
!_L    qn         : arg real(ip40      ) ; vitesse normale
!_L    p          : arg real(ip40      ) ; pression statique
!_L    resi       : arg real(ip40      ) ; residu ds resolution par Newton
!_L    ro         : arg real(ip11      ) ; masse volumique
!_L    un         : arg real(ip40      ) ; u.n, u vitesse d'entrainement(v abs=
!_L                                        w rel+u ent) et n normale a la front
!_L    usdn       : arg real(ip40      ) ; 1/(d.n), d direction de la vitesse
!_L                                        et n normale a la front
!_L    ym         : arg real(ip40      ) ; coordonnee sur l'axe y du centre
!_L                                        de la facette frontiere
!_L    zm         : arg real(ip40      ) ; coordonnee sur l'axe z du centre
!_L                                        de la facette frontiere
!
!     COM
!_C    Notation 0      : valeurs a l' instant n.
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
    integer          ::    id,  inc,   jd,   kd,    l
    integer          ::     m,   mb,  mfb,  mmb,   mn
    integer          ::   mpb,  mpn,   mt,  n0c,  n0n
    integer          ::    nc, ncbd,  nci, ncij,ncijk
    integer          ::  ncik,  ncj, ncjk,  nck,  nid
    integer          ::  nijd, nitn,  njd,   nn
    double precision ::       a2,     alm,     am0,       b,   beta2
    double precision ::     cson,     d0x,     d0y,     d0z,      df
    double precision ::      dpn,     eps,       f,   gam2t,    gam7
    double precision ::     gam8,     nxn,     nyn,     nzn,       p
    double precision :: pression,      ps,      q2,    qinf,      qn
    double precision ::      qxs,     qys,     qzs,    resi,  residu
    double precision ::       ro,    roc0,     rpi,     rti,    temp
    double precision ::       un,    usdn,   usdn2,       v,       w
    double precision ::       wn,       y,      ym,       z,      zm
!
!-----------------------------------------------------------------------
!
    dimension rpi(ip40),rti(ip40),d0x(ip40),d0y(ip40),d0z(ip40)
    dimension v(ip11,ip60)
    dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
    dimension mmb(mtt),mpb(mtt),mpn(mtt)
    dimension usdn2(ip40),roc0(ip40),am0(ip40),alm(ip40), &
         qn(ip40),p(ip40),resi(ip40),ro(ip40), &
         un(ip40),usdn(ip40),ym(ip40),zm(ip40)
    dimension pression(ip11),temp(ip11),cson(ip11)
    dimension y(ip21),z(ip21)
!
    inc(id,jd,kd)=id+jd*nid+kd*nijd
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
    gam7=-gam1/gam
    gam8=(1.-2.*gam)/gam
    gam2t=0.5*gam1/aa1**2
!
    qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
!
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
       usdn(m)=1./(d0x(m)*nxn(mn)+d0y(m)*nyn(mn)+d0z(m)*nzn(mn))
       usdn2(m)=usdn(m)**2
       un(m)=omg*(ym(m)*nzn(mn)-zm(m)*nyn(mn))
!
       qxs=v(nc,2)/v(nc,1)
       qys=v(nc,3)/v(nc,1)
       qzs=v(nc,4)/v(nc,1)
       qn(m)=qxs*nxn(mn)+qys*nyn(mn)+qzs*nzn(mn)
       q2=qxs**2+qys**2+qzs**2
       ps=gam1*(v(nc,5)-0.5*v(nc,1)*q2-pinfl)
       a2=gam*ps/v(nc,1)
       beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
       roc0(m)=v(nc,1)*beta2*a2
       alm(m)=0.5*((1.+beta2)*abs(qn(m))-sqrt(((1.-beta2)*qn(m))**2+4.*beta2*a2))
       am0(m)=(alm(m)-qn(m))*ps+roc0(m)*qn(m)
       p(m)=ps
    enddo
!
!   resolution de (alm-qn)*p + rho*beta2*a2*qn - am0 = 0 par Newton
!
    nitn=0
    residu=1.
    do m=1,mt
       mb=mpb(mfb)+m
       nc=ncbd(mb)
       nn=nc-n0c+n0n
       do while(residu.gt.eps)
          b=roc0(m)*sqrt(rti(m)/gam2t)/usdn(m)
          f=(alm(m)-qn(m))*p(m)+b*sqrt(abs((gam*p(m))**gam7-1.))-am0(m)
          df=alm(m)-qn(m) -0.5*gam1*b*((gam*p(m))**gam8)/ &
               sqrt(abs((gam*p(m))**gam7-1.))
          dpn=-f/df
          residu=abs(dpn)/p(m)
!         p(m)=p(m)+dpn
          p(m)=max(0.,p(m)+dpn)
          nitn=nitn+1
          if(nitn.gt.10) then
             residu=0.
             p(m)=gam1*(v(nc,5)-pinfl-0.5*(v(nc,2)**2+v(nc,3)**2+v(nc,4)**2)/v(nc,1))
          endif
       enddo
       nitn=0
       residu=1.
    enddo
!
    do m=1,mt
       mb=mpb(mfb)+m
       nc=ncbd(mb)
       qn(m)=sqrt(rti(m)/gam2t)/usdn(m)*sqrt((gam*p(m))**gam7 -1.)
       ro(m)=(gam*p(m))**(1./gam)
!       ro(m)=rti(m)
!       qn(m)=d0z(m)
       wn=qn(m)+un(m)
       w=wn*usdn(m)
!
       v(nc,1)=ro(m)
       v(nc,2)=ro(m)*(w*d0x(m))
       v(nc,3)=ro(m)*(w*d0y(m)+zm(m)*omg)
       v(nc,4)=ro(m)*(w*d0z(m)-ym(m)*omg)
       v(nc,5)=p(m)/gam1+pinfl +0.5*(v(nc,2)**2+v(nc,3)**2+v(nc,4)**2)/ro(m)
!
       pression(nc)=p(m)
       temp(nc)=gam*p(m)/ro(m)
       cson(nc)=sqrt(temp(nc))
    enddo
!
    return
  end subroutine clidd2_prcd
end module mod_clidd2_prcd
