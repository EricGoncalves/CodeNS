module mod_utitfr_xy
  implicit none
contains
  subroutine utitfr( &
       x,y,z,v,ncbd,nxn,nyn,nzn,icyc, &
       toxx,toxy,toxz,toyy,toyz,tozz, &
       ps)
!
!***********************************************************************
!
!     ACT
!_A    Sous-programme d'exploitation des resultats en cour de calcul.
!_A
!_A    calcul et ecriture dans fout des efforts de pression et frottement
!_A    de la configuration complete.
!_A    ATTENTION! demi-longueur en envergure.
!
!     VAL
!_V    Il faut que les surfaces projetees sur xy soient non nulles.
!
!     INP
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_I                                        contraintes visqueuses
!_I    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_I                                        contraintes visqueuses
!_I    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_I                                        contraintes visqueuses
!_I    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_I                                        contraintes visqueuses
!_I    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_I                                        contraintes visqueuses
!_I    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_I                                        contraintes visqueuses
!_I    it         : arg int              ; cycle courant du calcul
!_I    iminb      : arg int (mtt       ) ; indice min en i d'une frontiere
!_I    imaxb      : arg int (mtt       ) ; indice max en i d'une frontiere
!_I    jminb      : arg int (mtt       ) ; indice min en j d'une frontiere
!_I    jmaxb      : arg int (mtt       ) ; indice max en j d'une frontiere
!_I    kminb      : arg int (mtt       ) ; indice min en k d'une frontiere
!_I    kmaxb      : arg int (mtt       ) ; indice max en k d'une frontiere
!_I    out        : com int              ; unite logiq, moyennes des residus
!_I    npn        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab tous noeuds
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    ii2        : com int (lt        ) ; indice max en i reel
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jj2        : com int (lt        ) ; indice max en j reel
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    gam1       : com real             ; rap chal spec -1
!_I    pa1        : com real             ; pression d'arret de l'etat
!_I                                        de reference utilisateur adimensionne
!_I    xref       : com real             ; x du pt de ref pour calcul moments
!_I    yref       : com real             ; y du pt de ref pour calcul moments
!_I    zref       : com real             ; z du pt de ref pour calcul moments
!_I    sref       : com real             ; surface de ref pour calcul efforts et
!_I                                        moments
!_I    xlref      : com real             ; longueur ref pour calcul moments
!_I    alpha0     : com real             ; an:wgle d'incidence, def meca vol
!_I    beta0      : com real             ; angle de derapage, def meca vol
!_I    p0spi0     : com real             ; pres statique/pres d'arret, infini
!_I    q0spi0     : com real             ; pres dynamique {0.5*ro0*v0**2}/pres
!_I                                        d'arret, infini
!_I    kvglo      : com int              ; cle calcul des grandeurs globales
!_I    nbfll      : com int              ; nb de front d'integration des
!_I                                        grandeurs globales
!_I    nmfint     : com int (mtb       ) ; no des front d'integration des
!_I                                        grandeurs globales
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use definition
    use proprieteflu
    use sortiefichier
    use schemanum 
    use constantes
    implicit none
    integer          ::      i1,     i2,   icyc,icyexpl,   idf1
    integer          ::    idf2,  idfac,    idm,  imaxf,  iminf
    integer          ::      j1,     j2,  jmaxf,  jminf,     k1
    integer          ::      k2,  kmaxf,  kminf,      l,    m0b
    integer          ::     m0n,     m1,  m1max,m1maxm1,  m1min
    integer          ::      m2,  m2max,m2maxm1,  m2min,     mf
    integer          ::    mfac,  mfacn,    mfl,    n0c,    n0n
    integer          ::    ncbd,    nci,    ncj,    nck,  nfac1
    integer          ::   nfac2,  nfac3,  nfac4,  nfacf,    nid
    integer          ::    nijd,    njd
    double precision ::   alfar,  betar, claefr, claero,   clav
    double precision ::  clavfr,clavtfr,clavtot, cmaefr, cmaero
    double precision ::    cmav, cmavfr,cmavtfr,cmavtot, cnaefr
    double precision ::  cnaero,   cnav, cnavfr,cnavtfr,cnavtot
    double precision ::    csal,   csbe, cxaefr, cxaero,   cxav
    double precision ::  cxavfr,cxavtfr,cxavtot, cyaefr, cyaero
    double precision ::    cyav, cyavfr,cyavtfr,cyavtot, czaefr
    double precision ::  czaero,   czav, czavfr,czavtfr,czavtot
    double precision ::     dcl,  dclfr,    dcm,  dcmfr,    dcn
    double precision ::   dcnfr,    dcx,  dcxfr,    dcy,  dcyfr
    double precision ::     dcz,  dczfr,   dsml,   dsxy,   dsxz
    double precision ::    dsyz,    dx1,    dx2,    dy1,    dy2
    double precision ::     dz1,    dz2,    nxn,    nyn,    nzn
    double precision ::    pres,     ps,  pspi0,   snal,   snbe
    double precision :: taunorm,tauref0,   toxx,   toxy,   toxz
    double precision ::    toyy,   toyz,   tozz,    utx,   utxt
    double precision ::     uty,   utyt,    utz,   utzt,      v
    double precision ::       x,  xcfac,      y,  ycfac,      z
    double precision ::   zcfac
!
!-----------------------------------------------------------------------
!
!
    dimension x(ip21),y(ip21),z(ip21)
    dimension toxx(ip12),toxy(ip12),toxz(ip12),toyy(ip12),toyz(ip12),tozz(ip12)
    dimension v(ip11,ip60)
    dimension ncbd(ip41),ps(ip11)
    dimension nxn(ip42),nyn(ip42),nzn(ip42)
!
    icyexpl=mod(icyc,ncyexpl)
!
    if(kvglo.eq.0) return
    if(nbfll.eq.0) return
!
    alfar=alpha0/raddeg
    betar=beta0/raddeg
    csal=cos(alfar)
    snal=sin(alfar)
    csbe=cos(betar)
    snbe=sin(betar)
!
    if(abs(pa1).le.tiny(1.)) then
       write(imp,'("!!!utit: pa1=0, devient 1")')
       pa1=1.
    end if
    if(abs(q0spi0).le.tiny(1.)) then
       write(imp,'("!!!utit: q0spi0=0, devient 1")')
       q0spi0=1.
    endif
    if(abs(sref).le.tiny(1.)) then
       write(imp,'("!!!utit: sref=0, devient 1")')
       sref=1.
    endif
!
!     -------------------------------------------------------
!     adimensionnement du frottement
!
!     definitions dans "dfph.f" :
!     gam=GAMMA ; gam1=gam-1 ; gam2=.5*gam1
!
!     Coefficient de frottement parietal ramene aux grandeurs a l'infini amont:
!     Cf0=tau_p / 0.5 rho_0 U_0**2  : grandeurs avec dimensions
!
!     Le tenseur des contraintes toxx, etc... est ramene a:
!     rho_i  -> masse volumique conditions d'arret
!     a_i**2 -> vitesse du son dans les conditions d'arret au carre
!
!     Cf0   =   (tau_calcul) * tauref
!
!     tauref = 2 * rho_i * a_i**2 / ( rho_inf * V_inf**2 )
!                  <-----grandeurs avec dimension------>
!
    tauref0=2.*(1.+gam2*rm0**2)**(gam/gam1)/rm0**2
!
!     -------------------------------------------------------
!
    if(icyexpl.eq.0) then

!     pression
    cxavtot=0.
    cyavtot=0.
    czavtot=0.
    clavtot=0.
    cmavtot=0.
    cnavtot=0.
!     frottement
    cxavtfr=0.
    cyavtfr=0.
    czavtfr=0.
    clavtfr=0.
    cmavtfr=0.
    cnavtfr=0.
!
    do mf=1,nbfll
       mfl=nmfint(mf)
       l=ndlb(mfl)
!
       i1=ii1(l)
       i2=ii2(l)
       j1=jj1(l)
       j2=jj2(l)
       k1=kk1(l)
       k2=kk2(l)
       n0n=npn(l)
       n0c=npc(l)
!
       nid = id2(l)-id1(l)+1
       njd = jd2(l)-jd1(l)+1
       nijd = nid*njd
!
       nci=1
       ncj=nid
       nck=nijd
!
       m0b=mpb(mfl)
       m0n=mpn(mfl)
       iminf=iminb(mfl)
       imaxf=imaxb(mfl)
       jminf=jminb(mfl)
       jmaxf=jmaxb(mfl)
       kminf=kminb(mfl)
       kmaxf=kmaxb(mfl)
!
       m1min=1
       m2min=1
       if (iminf.eq.imaxf) then
          m1max=jmaxf-jminf+1
          m1maxm1=m1max-1
          m2max=kmaxf-kminf+1
          m2maxm1=m2max-1
          idf1=ncj
          idf2=nck
!
          if(iminf.eq.i1) idfac=nci
          if(iminf.eq.i2) idfac=0
!
       elseif (jminf.eq.jmaxf) then
          m1max=imaxf-iminf+1
          m1maxm1=m1max-1
          m2max=kmaxf-kminf+1
          m2maxm1=m2max-1
          idf1=nci
          idf2=nck
!
          if(jminf.eq.j1) idfac=ncj
          if(jminf.eq.j2) idfac=0
!
       elseif (kminf.eq.kmaxf) then
          m1max=imaxf-iminf+1
          m1maxm1=m1max-1
          m2max=jmaxf-jminf+1
          m2maxm1=m2max-1
          idf1=nci
          idf2=ncj
!
          if(kminf.eq.k1) idfac=nck
          if(kminf.eq.k2) idfac=0
!
       endif
!
       idm=m1max-m1min
!
!       pression
       cxav=0.
       cyav=0.
       czav=0.
       clav=0.
       cmav=0.
       cnav=0.
!       frottement
       cxavfr=0.
       cyavfr=0.
       czavfr=0.
       clavfr=0.
       cmavfr=0.
       cnavfr=0.
!
       do m2=m2min,m2maxm1
          do m1=m1min,m1maxm1
             mfac =m0b+m1+(m2-1)*idm
             mfacn=m0n+m1+(m2-1)*idm
             nfac1=ncbd(mfac)-n0c+n0n+idfac
             nfac2=ncbd(mfac)-n0c+n0n+idfac+idf1
             nfac3=ncbd(mfac)-n0c+n0n+idfac+idf1+idf2
             nfac4=ncbd(mfac)-n0c+n0n+idfac+idf2
             nfacf=ncbd(mfac)
!
!           pres=gam1*( v(nfacf,5) - pinfl &
!              -0.5*(v(nfacf,2)**2+v(nfacf,3)**2+v(nfacf,4)**2)/v(nfacf,1))
             pres=ps(nfacf)
             pspi0=pres/pa1          !pa1=1/gam
             xcfac=0.25*(x(nfac1)+x(nfac2)+x(nfac3)+x(nfac4))
             ycfac=0.25*(y(nfac1)+y(nfac2)+y(nfac3)+y(nfac4))
             zcfac=0.25*(z(nfac1)+z(nfac2)+z(nfac3)+z(nfac4))
             dx1=x(nfac3)-x(nfac1)
             dy1=y(nfac3)-y(nfac1)
             dz1=z(nfac3)-z(nfac1)
             dx2=x(nfac4)-x(nfac2)
             dy2=y(nfac4)-y(nfac2)
             dz2=z(nfac4)-z(nfac2)
             dsyz=abs(dy1*dz2-dz1*dy2)/2.
             dsxz=abs(dz1*dx2-dx1*dz2)/2.
             dsxy=abs(dx1*dy2-dy1*dx2)/2.
             dsml=sqrt(dsyz*dsyz+dsxz*dsxz+dsxy*dsxy)
!           pression
             dcx=(p0spi0-pspi0)*dsyz*sign(1.D0,nxn(mfacn))
             dcy=(p0spi0-pspi0)*dsxz*sign(1.D0,nyn(mfacn))
             dcz=(p0spi0-pspi0)*dsxy*sign(1.D0,nzn(mfacn))
             dcl=dcy*(zcfac-zref)-dcz*(ycfac-yref)
             dcm=dcx*(zcfac-zref)-dcz*(xcfac-xref)
             dcn=dcx*(ycfac-yref)-dcy*(xcfac-xref)
             cxav=cxav+dcx
             cyav=cyav+dcy
             czav=czav+dcz
             clav=clav+dcl
             cmav=cmav+dcm
             cnav=cnav+dcn
!
!           frottement
             utx=toxx(nfacf)*nxn(mfacn)+toxy(nfacf)*nyn(mfacn)+ &
                 toxz(nfacf)*nzn(mfacn)
             uty=toxy(nfacf)*nxn(mfacn)+toyy(nfacf)*nyn(mfacn)+ &
                 toyz(nfacf)*nzn(mfacn)
             utz=toxz(nfacf)*nxn(mfacn)+toyz(nfacf)*nyn(mfacn)+ &
                 tozz(nfacf)*nzn(mfacn)
!           projection du frottement sur la surface
             taunorm=utx*nxn(mfacn)+uty*nyn(mfacn)+utz*nzn(mfacn)
             utxt=utx-taunorm*nxn(mfacn)
             utyt=uty-taunorm*nyn(mfacn)
             utzt=utz-taunorm*nzn(mfacn)
             dcxfr=utxt*dsml
             dcyfr=utyt*dsml
             dczfr=utzt*dsml
             dclfr=dcyfr*(zcfac-zref)-dczfr*(ycfac-yref)
             dcmfr=dcxfr*(zcfac-zref)-dczfr*(xcfac-xref)
             dcnfr=dcxfr*(ycfac-yref)-dcyfr*(xcfac-xref)
             cxavfr=cxavfr+dcxfr
             cyavfr=cyavfr+dcyfr
             czavfr=czavfr+dczfr
             clavfr=clavfr+dclfr
             cmavfr=cmavfr+dcmfr
             cnavfr=cnavfr+dcnfr
          enddo
       enddo
!
!      pression
       cxav=cxav/(q0spi0*sref)
       cyav=cyav/(q0spi0*sref)
       czav=czav/(q0spi0*sref)
       clav=clav/(q0spi0*sref*xlref)
       cmav=cmav/(q0spi0*sref*xlref)
       cnav=cnav/(q0spi0*sref*xlref)
!       cxaero=cxav*csal*csbe-cyav*snbe+czav*snal*csbe
!       cyaero=cxav*csal*snbe+cyav*csbe+czav*snal*snbe
!       czaero=-cxav*snal+czav*csal
       cxaero=cxav*csal*csbe-czav*snbe+cyav*snal*csbe
       czaero=cxav*csal*snbe+czav*csbe+cyav*snal*snbe
       cyaero=-cxav*snal+cyav*csal
       claero=clav*csal*csbe+cmav*snbe+cnav*snal*csbe
       cmaero=-clav*csal*snbe+cmav*csbe-cnav*snal*snbe
       cnaero=-clav*snal+cnav*csal
       cxavtot=cxavtot+cxav
       cyavtot=cyavtot+cyav
       czavtot=czavtot+czav
       clavtot=clavtot+clav
       cmavtot=cmavtot+cmav
       cnavtot=cnavtot+cnav
!
!      frottement
       cxavfr=cxavfr*tauref0/sref
       cyavfr=cyavfr*tauref0/sref
       czavfr=czavfr*tauref0/sref
       clavfr=clavfr*tauref0/(sref*xlref)
       cmavfr=cmavfr*tauref0/(sref*xlref)
       cnavfr=cnavfr*tauref0/(sref*xlref)
!
       cxavtfr=cxavtfr+cxavfr
       cyavtfr=cyavtfr+cyavfr
       czavtfr=czavtfr+czavfr
       clavtfr=clavtfr+clavfr
       cmavtfr=cmavtfr+cmavfr
       cnavtfr=cnavtfr+cnavfr
    enddo !fin boucle sur les parois
!
!   pression
!    cxaero= cxavtot*csal*csbe-cyavtot*snbe+czavtot*snal*csbe
!    cyaero= cxavtot*csal*snbe+cyavtot*csbe+czavtot*snal*snbe
!    czaero=-cxavtot*snal+czavtot*csal
    cxaero= cxavtot*csal*csbe-czavtot*snbe+cyavtot*snal*csbe
    czaero= cxavtot*csal*snbe+czavtot*csbe+cyavtot*snal*snbe
    cyaero=-cxavtot*snal+cyavtot*csal
    claero= clavtot*csal*csbe+cmavtot*snbe+cnavtot*snal*csbe
    cmaero=-clavtot*csal*snbe+cmavtot*csbe-cnavtot*snal*snbe
    cnaero=-clavtot*snal+cnavtot*csal
!
!   frottement
!    cxaefr= cxavtfr*csal*csbe-cyavtfr*snbe+czavtfr*snal*csbe
!    cyaefr= cxavtfr*csal*snbe+cyavtfr*csbe+czavtfr*snal*snbe
!    czaefr=-cxavtfr*snal+czavtfr*csal
    cxaefr= cxavtfr*csal*csbe-czavtfr*snbe+cyavtfr*snal*csbe
    czaefr= cxavtfr*csal*snbe+czavtfr*csbe+cyavtfr*snal*snbe
    cyaefr=-cxavtfr*snal+cyavtfr*csal
    claefr= clavtfr*csal*csbe+cmavtfr*snbe+cnavtfr*snal*csbe
    cmaefr=-clavtfr*csal*snbe+cmavtfr*csbe-cnavtfr*snal*snbe
    cnaefr=-clavtfr*snal+cnavtfr*csal

!    if(icyexpl.eq.0) then
!      repere avion
       write(out,3801) icyc,cxavtot,cyavtot,czavtot, &
            clavtot,cmavtot,cnavtot
3801   format('=>utitfr p: ',i6,1x,6(1pe11.3))
       write(out,3802) icyc,cxavtfr,cyavtfr,czavtfr,clavtfr, &
            cmavtfr,cnavtfr
3802   format('=>utitfr f: ',i6,1x,6(1pe11.3))
    endif
!
!   blocage de la condition vrt
!   commente la ligne pour bloquer 
!    vrtcz=czaero
!
    return
  end subroutine utitfr
end module mod_utitfr_xy
