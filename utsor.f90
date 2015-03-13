module mod_utsor
implicit none
contains
      subroutine utsor( &
                 ncbd,s,mut, &
                 x,y,z,nxn,nyn,nzn)
!
!***********************************************************************
!
!     ACT
!_A    Sous-programme utilisateur d'exploitation des resultats a la fin
!_A    du calcul.
!_A
!_A     Ecriture dans fsor1 des informations necessaires a
!_A     l'exploitation de resultats ( fichiers fg + fa ) a
!_A     l'exterieur du programme principal :
!_A        -geometrie
!_A        -physique
!_A        -ecoulement
!_A        -normalisation
!_A        -etat de reference
!_A     Dans le cas ou la cle 'kvglo' est non nulle,
!_A     calcul et ecriture dans fimp de
!_A        -description des surfaces d'integration
!_A        -efforts par bandes
!_A        -efforts par surfaces
!_A        -efforts de la configuration complete
!_A     calcul et ecriture dans fsor2 pour chaque facette
!_A        -des coordonnees du centre de la facette
!_A        -des composantes de la vitesse rapportees a v0
!_A        -du kp {(p-p0)/(.5*ro0*v0**2)}
!_A        -de pi/pi0
!
!     VAL
!_V    Il faut que les surfaces projetees sur xy
!_V    soient non nulles.
!
!     INP
!_I    ip11       : arg int              ; dim, nbr max de cellules de tous les
!_I                                        dom (pts fictifs inclus)
!_I    ip60       : arg int              ; dim, nbr max d'equations
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    s          : arg real(ip11,ip60 ) ; variables de calcul
!_I    omg        : arg real             ; vitesse rotation du repere relatif
!_I    perio      : arg real             ; periodicite geometrique en angle ou
!_I                                        distance selon config
!_I    kimp       : arg int              ; niveau de sortie sur unite logi imp
!_I    equat      : arg char             ; type d'equations modelisant l'ecoule-
!_I                                        ment
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    sorf1      : com int              ; unite logiq, sorties utilisateur 1
!_I    sorf2      : com int              ; unite logiq, sorties utilisateur 2
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
!_I    iminb      : com int (mtt       ) ; indice min en i d'une frontiere
!_I    imaxb      : com int (mtt       ) ; indice max en i d'une frontiere
!_I    jminb      : com int (mtt       ) ; indice min en j d'une frontiere
!_I    jmaxb      : com int (mtt       ) ; indice max en j d'une frontiere
!_I    kminb      : com int (mtt       ) ; indice min en k d'une frontiere
!_I    kmaxb      : com int (mtt       ) ; indice max en k d'une frontiere
!_I    roa1       : com real             ; etat de reference utilisateur
!_I                                        adimensionne, masse volumique d'arret
!_I    aa1        : com real             ; etat de reference utilisateur
!_I                                        adimensionne, vitesse du son d'arret
!_I    ta1        : com real             ; etat de reference utilisateur
!_I                                        adimensionne, temperature d'arret
!_I    pa1        : com real             ; pression d'arret de l'etat
!_I                                        de reference utilisateur adimensionne
!_I    ha1        : com real             ; enthalpie d'arret de l'etat
!_I                                        de reference utilisateur adimensionne
!_I    gam        : com real             ; rapport des chaleurs specifiques
!_I    gam1       : com real             ; rap chal spec -1
!_I    gam2       : com real             ; (rap chal spec -1)/2
!_I    gam4       : com real             ; 1/(rap chal spec -1)
!_I    cp         : com real             ; chal spec a pres cste adim
!_I    cv         : com real             ; chal spec a vol cst adim
!_I    pr         : com real             ; nombre de Prandtl
!_I    prt        : com real             ; nombre de Prandtl turbulent
!_I    reynz      : com real             ; nombre de Reynolds calcule avec
!_I                                        les grandeurs d'adimensionnement,
!_I                                        pour definir la loi de Sutherland
!_I    tnz        : com real             ; etat pour adimensionnement,
!_I                                        temperature
!_I    ronz       : com real             ; etat pour adimensionnement,
!_I                                        masse volumique
!_I    anz        : com real             ; etat pour adimensionnement,
!_I                                        vitesse du son d'arret
!_I    dnz        : com real             ; etat pour adimensionnement,
!_I                                        longueur
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    xref       : com real             ; x du pt de ref pour calcul grandeurs
!_I                                        globales
!_I    yref       : com real             ; y du pt de ref pour calcul grandeurs
!_I                                        globales
!_I    zref       : com real             ; z du pt de ref pour calcul grandeurs
!_I                                        globales
!_I    sref       : com real             ; surface de ref pour calcul grandeurs
!_I                                        globales
!_I    xlref      : com real             ; longueur de ref pour calcul grandeurs
!_I                                        globales
!_I    alpha0     : com real             ; angle d'incidence
!_I    beta0      : com real             ; angle de derapage
!_I    p0spi0     : com real             ; pres statique/pres d'arret, infini
!_I    q0spi0     : com real             ; vitesse/pres d'arret, infini
!_I    v0         : com real             ; vitesse de l'ecoulement a l'infini
!_I    kvglo      : com int              ; cle calcul dea grandeurs globales
!_I    nbfll      : com int              ; nb de front d'integration des
!_I                                        grandeurs globales
!_I    nmfint     : com int (mtb       ) ; no des front d'integration des
!_I                                        grandeurs globales
!
!
!     COM
!_C    Les efforts et les moments sont calcules a la fois dans le repere
!_C    avion et dans le repere aerodynamique
!
!***********************************************************************
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use boundary
      use sortiefichier
      use proprieteflu 
      use definition
      use chainecarac
implicit none
integer :: ncbd
double precision :: s
double precision :: x
double precision :: y
double precision :: z
double precision :: akp
double precision :: alfar
double precision :: betar
double precision :: claero
double precision :: claerob
double precision :: clav
double precision :: clavb
double precision :: clavtot
double precision :: cmaero
double precision :: cmaerob
double precision :: cmav
double precision :: cmavb
double precision :: cmavtot
double precision :: cnaero
double precision :: cnaerob
double precision :: cnav
double precision :: cnavb
double precision :: cnavtot
double precision :: csal
double precision :: csbe
double precision :: cxaero
double precision :: cxaerob
double precision :: cxav
double precision :: cxavb
double precision :: cxavtot
double precision :: cyaero
double precision :: cyaerob
double precision :: cyav
double precision :: cyavb
double precision :: cyavtot
double precision :: czaero
double precision :: czaerob
double precision :: czav
double precision :: czavb
double precision :: czavtot
double precision :: dcl
double precision :: dcm
double precision :: dcn
double precision :: dcx
double precision :: dcy
double precision :: dcz
double precision :: dsml
double precision :: dsxy
double precision :: dsxz
double precision :: dsyz
double precision :: dx1
double precision :: dx2
double precision :: dy1
double precision :: dy2
double precision :: dz1
double precision :: dz2
integer :: i1
integer :: i2
integer :: idf1
integer :: idf2
integer :: idfac
integer :: idm
integer :: imaxf
integer :: iminf
integer :: j1
integer :: j2
integer :: jmaxf
integer :: jminf
integer :: k1
integer :: k2
integer :: kmaxf
integer :: kminf
integer :: l
integer :: m0b
integer :: m0n
integer :: m1
integer :: m1max
integer :: m1maxm1
integer :: m1min
integer :: m2
integer :: m2max
integer :: m2maxm1
integer :: m2min
integer :: mf
integer :: mfac
integer :: mfacn
integer :: mfl
integer :: n0c
integer :: n0n
integer :: nci
integer :: ncj
integer :: nck
integer :: nfac1
integer :: nfac2
integer :: nfac3
integer :: nfac4
integer :: nfacf
integer :: nid
integer :: nijd
integer :: njd
double precision :: pa
double precision :: paspa1
double precision :: pis2
double precision :: ps
double precision :: pspi0
double precision :: qq
double precision :: raddeg
double precision :: rm2
double precision :: sml
double precision :: snal
double precision :: snbe
double precision :: sxy
double precision :: sxyb
double precision :: sxz
double precision :: syz
double precision :: u
double precision :: v
double precision :: w
double precision :: xcfac
double precision :: ycfac
double precision :: zcfac
!
!-----------------------------------------------------------------------
!
      double precision mut,nxn,nyn,nzn
!
      dimension x(ip21),y(ip21),z(ip21)
      dimension s(ip11,ip60)
      dimension mut(ip12),nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
!
!     SORTIES POUR EXPLOITATION
!
      open(sorf1,file='fsor1')
!
      write(sorf1,1000) equat
!
      if(equat(1:2).eq.'eu') then
      write(sorf1,2000) cp,cv
      elseif(equat(1:2).eq.'ns') then
      write(sorf1,2000) cp,cv,pr,prt,reynz
      endif
!
      write(sorf1,2000) tnz,ronz,anz,dnz
      write(sorf1,2000) roa1,aa1,ta1,pa1,ha1
      write(sorf1,2000) omg,perio
!
 1000 format(a)
 2000 format(5e14.7)
!
!     SORTIES RELATIVES A DES VALEURS SUR LES PAROIS
!
      if(kvglo.eq.0) return
      if(nbfll.eq.0) return
!
      open(sorf2,file='fsor2')
!
      pis2=atan2(1.,0.)
      raddeg=90./pis2
!
      alfar=alpha0/raddeg
      betar=beta0/raddeg
      csal=cos(alfar)
      snal=sin(alfar)
      csbe=cos(betar)
      snbe=sin(betar)
!
      cxavtot=0.
      cyavtot=0.
      czavtot=0.
      clavtot=0.
      cmavtot=0.
      cnavtot=0.
!
      do mf=1,nbfll
!
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
      if(kimp.eq.1) then
       write(imp,987) l,mfl,iminf,imaxf,jminf,jmaxf,kminf,kmaxf
      endif
  987 format("integration des pressions par frontiere : ",39("-") &
      //1x,"zone ",i3,"   - frontiere ",i3/1x,26("-"), &
      //5x,"imin = ",i3,5x,"imax = ",i3, &
      //5x,"jmin = ",i3,5x,"jmax = ",i3, &
      //5x,"kmin = ",i3,5x,"kmax = ",i3/)
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
      cxav=0.
      cyav=0.
      czav=0.
      clav=0.
      cmav=0.
      cnav=0.
      syz=0.
      sxz=0.
      sxy=0.
      sml=0.
!
      do m2=m2min,m2maxm1
      cxavb=0.
      cyavb=0.
      czavb=0.
      clavb=0.
      cmavb=0.
      cnavb=0.
      sxyb=0.
!
      do m1=m1min,m1maxm1
      mfac=m0b+m1+(m2-1)*idm
      mfacn=m0n+m1+(m2-1)*idm
      nfac1=ncbd(mfac)-n0c+n0n+idfac
      nfac2=ncbd(mfac)-n0c+n0n+idfac+idf1
      nfac3=ncbd(mfac)-n0c+n0n+idfac+idf1+idf2
      nfac4=ncbd(mfac)-n0c+n0n+idfac+idf2
      nfacf=ncbd(mfac)
      ps =gam1*( s(nfacf,1)*s(nfacf,5) &
           -0.5*(s(nfacf,2)**2+s(nfacf,3)**2+s(nfacf,4)**2))/s(nfacf,1)
      pspi0=ps/pa1
      xcfac=(x(nfac1)+x(nfac2)+x(nfac3)+x(nfac4))/4.
      ycfac=(y(nfac1)+y(nfac2)+y(nfac3)+y(nfac4))/4.
      zcfac=(z(nfac1)+z(nfac2)+z(nfac3)+z(nfac4))/4.
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
      dcx=(p0spi0-pspi0)*dsyz*sign(1.,nxn(mfacn))
      dcy=(p0spi0-pspi0)*dsxz*sign(1.,nyn(mfacn))
      dcz=(p0spi0-pspi0)*dsxy*sign(1.,nzn(mfacn))
      dcl=dcy*(zcfac-zref)-dcz*(ycfac-yref)
      dcm=dcx*(zcfac-zref)-dcz*(xcfac-xref)
      dcn=dcx*(ycfac-yref)-dcy*(xcfac-xref)
      cxav=cxav+dcx
      cyav=cyav+dcy
      czav=czav+dcz
      clav=clav+dcl
      cmav=cmav+dcm
      cnav=cnav+dcn
      syz=syz+dsyz
      sxz=sxz+dsxz
      sxy=sxy+dsxy
      sml=sml+dsml
      cxavb=cxavb+dcx
      cyavb=cyavb+dcy
      czavb=czavb+dcz
      clavb=clavb+dcl
      cmavb=cmavb+dcm
      cnavb=cnavb+dcn
      sxyb=sxyb+dsxy
      enddo
!
      if(sxyb.gt.0.) then
      cxavb=cxavb/(q0spi0*sxyb)
      cyavb=cyavb/(q0spi0*sxyb)
      czavb=czavb/(q0spi0*sxyb)
      clavb=clavb/(q0spi0*sxyb*xlref)
      cmavb=cmavb/(q0spi0*sxyb*xlref)
      cnavb=cnavb/(q0spi0*sxyb*xlref)
      cxaerob=cxavb*csal*csbe-cyavb*snbe+czavb*snal*csbe
      cyaerob=cxavb*csal*snbe+cyavb*csbe+czavb*snal*snbe
      czaerob=-cxavb*snal+czavb*csal
      claerob=clavb*csal*csbe+cmavb*snbe+cnavb*snal*csbe
      cmaerob=-clavb*csal*snbe+cmavb*csbe-cnavb*snal*snbe
      cnaerob=-clavb*snal+cnavb*csal
!
      if(kimp.eq.1) then
      write(imp,991) m2,sxyb,cxavb,cyavb,czavb,clavb,cmavb,cnavb, &
      cxaerob,cyaerob,czaerob,claerob,cmaerob,cnaerob
      endif
  991 format(//,1x,"bande numero ",i3," :",/,1x,18("-"),// &
      /,1x,"surface mouillee projetee sur xy : ",e12.4,/ &
      /,1x,"efforts dans le repere avion : ",/ &
      /,5x,"cx = ",f8.4,5x,"cy = ",f8.4,5x,"cz = ",f8.4 &
      ,5x,"cl = ",f8.4,5x,"cm = ",f8.4,5x,"cn = ",f8.4,// &
      /,1x,"efforts dans le repere aerodynamique : ",/ &
      /,5x,"cx = ",f8.4,5x,"cy = ",f8.4,5x,"cz = ",f8.4 &
      ,5x,"cl = ",f8.4,5x,"cm = ",f8.4,5x,"cn = ",f8.4//)
!
      endif
!
!  ecriture des variables aux centres des mailles surfaciques :
!
!     x, y, z, u/v0, v/v0, w/v0, kp, pa/pa1
!
      do m1=m1min,m1maxm1
      mfac=m0b+m1+(m2-1)*idm
      nfac1=ncbd(mfac)-n0c+n0n+idfac
      nfac2=ncbd(mfac)-n0c+n0n+idfac+idf1
      nfac3=ncbd(mfac)-n0c+n0n+idfac+idf1+idf2
      nfac4=ncbd(mfac)-n0c+n0n+idfac+idf2
      nfacf=ncbd(mfac)
      ps =gam1*( s(nfacf,1)*s(nfacf,5) &
           -0.5*(s(nfacf,2)**2+s(nfacf,3)**2+s(nfacf,4)**2)) &
               /s(nfacf,1)
      pspi0=ps/pa1
      xcfac=(x(nfac1)+x(nfac2)+x(nfac3)+x(nfac4))/4.
      ycfac=(y(nfac1)+y(nfac2)+y(nfac3)+y(nfac4))/4.
      zcfac=(z(nfac1)+z(nfac2)+z(nfac3)+z(nfac4))/4.
      u=s(nfacf,2)/(s(nfacf,1)*v0)
      v=s(nfacf,3)/(s(nfacf,1)*v0)
      w=s(nfacf,4)/(s(nfacf,1)*v0)
      qq=(u*u+v*v+w*w)*v0*v0
      rm2=qq*s(nfacf,1)/(gam*ps)
      pa=ps*(1.+gam2*rm2)**(gam*gam4)
      paspa1=pa/pa1
      akp=(pspi0-p0spi0)/q0spi0
      write(sorf2,555) xcfac,ycfac,zcfac,u,v,w,akp,paspa1
 555  format(8f10.6)
      enddo
!
      enddo
!
      cxav=cxav/(q0spi0*sref)
      cyav=cyav/(q0spi0*sref)
      czav=czav/(q0spi0*sref)
      clav=clav/(q0spi0*sref*xlref)
      cmav=cmav/(q0spi0*sref*xlref)
      cnav=cnav/(q0spi0*sref*xlref)
      cxaero=cxav*csal*csbe-cyav*snbe+czav*snal*csbe
      cyaero=cxav*csal*snbe+cyav*csbe+czav*snal*snbe
      czaero=-cxav*snal+czav*csal
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
      if(kimp.eq.1) then
       write(imp,985)
       write(imp,989) sml,sxy,syz,sxz,cxav,cyav,czav,clav,cmav, &
       cnav,cxaero,cyaero,czaero,claero,cmaero,cnaero
      endif
  985 format(//,1x,"frontiere complete :",/,1x,20("-"),/)
  989 format(/,1x,"surface mouillee          : ",e12.4,/ &
      /,1x,"surface mouillee projetee sur xy : ",e12.4,/ &
      /,1x,"surface mouillee projetee sur yz : ",e12.4,/ &
      /,1x,"surface mouillee projetee sur xz : ",e12.4,// &
      /,1x,"efforts dans le repere avion : ",/ &
      /,5x,"cx = ",f8.4,5x,"cy = ",f8.4,5x,"cz = ",f8.4 &
      ,5x,"cl = ",f8.4,5x,"cm = ",f8.4,5x,"cn = ",f8.4,// &
      /,1x,"efforts dans le repere aerodynamique : ",/ &
      /,5x,"cx = ",f8.4,5x,"cy = ",f8.4,5x,"cz = ",f8.4 &
      ,5x,"cl = ",f8.4,5x,"cm = ",f8.4,5x,"cn = ",f8.4//)
      enddo
!
      cxaero=cxavtot*csal*csbe-cyavtot*snbe+czavtot*snal*csbe
      cyaero=cxavtot*csal*snbe+cyavtot*csbe+czavtot*snal*snbe
      czaero=-cxavtot*snal+czavtot*csal
      claero=clavtot*csal*csbe+cmavtot*snbe+cnavtot*snal*csbe
      cmaero=-clavtot*csal*snbe+cmavtot*csbe-cnavtot*snal*snbe
      cnaero=-clavtot*snal+cnavtot*csal
!
      if(kimp.eq.1) then
      write(imp,990) cxavtot,cyavtot,czavtot,clavtot,cmavtot,cnavtot, &
      cxaero,cyaero,czaero,claero,cmaero,cnaero
      endif
  990 format("efforts globaux - configuration complete :"/,1x, &
      40("-"),//,1x,"repere avion :",/ &
      /5x,"cx = ",f8.4,5x,"cy = ",f8.4,5x,"cz = ",f8.4 &
      ,5x,"cl = ",f8.4,5x,"cm = ",f8.4,5x,"cn = ",f8.4,// &
      /,1x,"repere aerodynamique :",/ &
      /,5x,"cx = ",f8.4,5x,"cy = ",f8.4,5x,"cz = ",f8.4 &
      ,5x,"cl = ",f8.4,5x,"cm = ",f8.4,5x,"cn = ",f8.4)
!
      return
      end subroutine
end module
