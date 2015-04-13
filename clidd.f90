module mod_clidd
  implicit none
contains
  subroutine clidd( &
       mfb,l,rpi,rti,d0x,d0y,d0z, &
       nxn,nyn,nzn,ncbd,v, &
       y,z, &
       mmb,mpb,mpn, &
       usdn2,roc0,am0,qn,p,resi,ro,un,usdn, &
       ym,zm,pression,temp,cson)
!
!***********************************************************************
!
!_DA  DATE_C : avril 2002  - Eric GONCALVES / SINUMEF
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement des conditions aux limites :
!_A     direction de la vitesse absolue imposee, parallele a (d0x,d0y,d0z),
!_A     pression d'arret imposee , pi=pa1*rpi
!_A     enthalpie d'arret imposee ,  hi=ha1*rti
!_A    Utilisation d une relation de compatibilite discretisee associee
!_A    a la caracteristique sortante (lambda=u-c)
!_A     Dp - rho*c*Dun = 0
!_A    Normales interieures.
!
!     INP
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
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    omg        : arg real             ; vitesse rotation du repere relatif
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    gam        : com real             ; rapport des chaleurs specifiques
!_I    gam1       : com real             ; rap chal spec -1
!_I    gam2       : com real             ; (rap chal spec -1)/2
!_I    gam4       : com real             ; 1/(rap chal spec -1)
!_I    gam5       : com real             ; rap chal spec/(rap chal spec -1)
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
    use definition
    implicit none
    integer          ::         id,        jd,        kd,         l,         m
    integer          ::         mb,       mfb,  mmb(mtt),        mn,  mpb(mtt)
    integer          ::   mpn(mtt),        mt,       n0c,       n0n,        nc
    integer          :: ncbd(ip41),       nci,      ncij,     ncijk,      ncik
    integer          ::        ncj,      ncjk,       nck,       nid,      nijd
    integer          ::       nitn,       njd,        nn
    double precision ::      am0(ip40),             b,    cson(ip11),     d0x(ip40),     d0y(ip40)
    double precision ::      d0z(ip40),            df,           dqn,           eps,             f
    double precision ::          gam2t,          gam6,          gamt,     nxn(ip42),     nyn(ip42)
    double precision ::      nzn(ip42),       p(ip40),pression(ip11),            ps,      qn(ip40)
    double precision ::            qxs,           qys,           qzs,    resi(ip40),        residu
    double precision ::       ro(ip40),    roc0(ip40),     rpi(ip40),     rti(ip40),    temp(ip11)
    double precision ::       un(ip40),    usdn(ip40),   usdn2(ip40),  v(ip11,ip60),             w
    double precision ::             wn,       y(ip21),      ym(ip40),       z(ip21),      zm(ip40)
!$OMP MASTER
!
!-----------------------------------------------------------------------
!
!
!

!
    eps=0.0000001
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
    gam2t=0.5*gam1/aa1**2
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
       usdn(m)=1./(d0x(m)*nxn(mn)+d0y(m)*nyn(mn)+d0z(m)*nzn(mn))
       usdn2(m)=usdn(m)**2
       un(m)=omg*(ym(m)*nzn(mn)-zm(m)*nyn(mn))
!
       qxs=v(nc,2)/v(nc,1)
       qys=v(nc,3)/v(nc,1)
       qzs=v(nc,4)/v(nc,1)
       qn(m)=qxs*nxn(mn)+qys*nyn(mn)+qzs*nzn(mn)
       ps=gam1*(v(nc,5)-0.5*(qxs*v(nc,2)+qys*v(nc,3)+qzs*v(nc,4))-pinfl)
       roc0(m)=-sqrt(gam*ps*v(nc,1))                  !attention au sign
       am0(m)=ps+roc0(m)*qn(m)
    enddo
!
!   resolution de p - rho*c*qn -am0 = 0 par algorithme de Newton
!
    nitn=0
    residu=1.
    do while( (residu.gt.eps).and.(nitn.lt.20) )
       residu=0.
!$OMP SIMD
       do m=1,mt
          mb=mpb(mfb)+m
!------la notation b designe la quantite 1/(1+((gam-1)/2)m**2)
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
       wn=qn(m)+un(m)
       w=wn*usdn(m)
!
       v(nc,1)=ro(m)
       v(nc,2)=ro(m)*(w*d0x(m))
       v(nc,3)=ro(m)*(w*d0y(m)+zm(m)*omg)
       v(nc,4)=ro(m)*(w*d0z(m)-ym(m)*omg)
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
  end subroutine clidd
end module mod_clidd
