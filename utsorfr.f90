module mod_utsorfr
  implicit none
contains
  subroutine utsorfr( &
       ncbd,ncin,s,mu,mut, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       x,y,z,nxn,nyn,nzn, &
       ps,temp,cson)
!
!***********************************************************************
!
!     ACT
!_A    Sous-programme d'exploitation des resultats a la fin du calcul.
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
!_A     calcul et ecriture dans pres pour chaque facette (frt)
!_A        -des coordonnees du centre de la facette
!_A        -des composantes du frottement parietal tau_p/(rho0 U0**2)
!_A        -du coefficient  de frottement parietal tau_p/(rho0 U0**2)
!_A        -du kp {(p-p0)/(.5*rho0*v0**2)}
!_A        -de pi/pi0
!
!     VAL
!_V    Il faut que les surfaces projetees sur xy soient non nulles.
!_A    ATTENTION! demi-longueur en envergure.
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
!_I    qcx        : arg real(ip12      ) ; composante en x du flux de chaleur
!_I    qcy        : arg real(ip12      ) ; composante en y du flux de chaleur
!_I    qcz        : arg real(ip12      ) ; composante en z du flux de chaleur
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
    use proprieteflu
    use definition
    use chainecarac
    use sortiefichier
    use schemanum
    use mod_mpi
    implicit none
    integer          ::         i1,        i2,      idf1,      idf2,     idfac
    integer          ::        idm,     imaxf,     iminf,        j1,        j2
    integer          ::      jmaxf,     jminf,        k1,        k2,     kimpl
    integer          ::      kmaxf,     kminf,         l,       m0b,       m0n
    integer          ::         m1,     m1max,   m1maxm1,     m1min,        m2
    integer          ::      m2max,   m2maxm1,     m2min,        mf,      mfac
    integer          ::      mfacn,       mfl,       n0c,       n0n,ncbd(ip41)
    integer          ::        nci,ncin(ip41),       ncj,       nck,     nfac1
    integer          ::      nfac2,     nfac3,     nfac4,     nfacf,     nfaci
    integer          ::      nfacm,       nid,      nijd,       njd,        nn,mfg
    double precision ::          akp,       alfar,       betar,      cfinf0,      claefr
    double precision ::      claerfr,      claero,     claerob,     claetfr,        clav
    double precision ::        clavb,     clavbfr,      clavfr,     clavtfr,     clavtot
    double precision ::       cmaefr,     cmaerfr,      cmaero,     cmaerob,     cmaetfr
    double precision ::         cmav,       cmavb,     cmavbfr,      cmavfr,     cmavtfr
    double precision ::      cmavtot,      cnaefr,     cnaerfr,      cnaero,     cnaerob
    double precision ::      cnaetfr,        cnav,       cnavb,     cnavbfr,      cnavfr
    double precision ::      cnavtfr,     cnavtot,        csal,        csbe,  cson(ip11)
    double precision ::       cxaefr,     cxaerfr,      cxaero,     cxaerob,     cxaetfr
    double precision ::         cxav,       cxavb,     cxavbfr,      cxavfr,     cxavtfr
    double precision ::      cxavtot,      cyaefr,     cyaerfr,      cyaero,     cyaerob
    double precision ::      cyaetfr,        cyav,       cyavb,     cyavbfr,      cyavfr
    double precision ::      cyavtfr,     cyavtot,      czaefr,     czaerfr,      czaero
    double precision ::      czaerob,     czaetfr,        czav,       czavb,     czavbfr
    double precision ::       czavfr,     czavtfr,     czavtot,         dcl,       dclfr
    double precision ::          dcm,       dcmfr,         dcn,       dcnfr,         dcx
    double precision ::        dcxfr,         dcy,       dcyfr,         dcz,       dczfr
    double precision ::         dsml,        dsxy,        dsxz,        dsyz,         dx1
    double precision ::          dx2,         dy1,         dy2,         dz1,         dz2
    double precision ::     mu(ip12),   mut(ip12),   nxn(ip42),   nyn(ip42),   nzn(ip42)
    double precision ::         phip,     phipadi,     phiref0,        phit,        pis2
    double precision ::     ps(ip11),       pspi0,   qcx(ip12),   qcy(ip12),   qcz(ip12)
    double precision ::           qq,          qt,         qtx,         qty,         qtz
    double precision ::       raddeg,          rm,s(ip11,ip60),         sml,       smlfr
    double precision ::         snal,        snbe,         sxy,        sxyb,       sxyfr
    double precision ::          sxz,       sxzfr,         syz,       syzfr,          t1
    double precision ::           t2,          t3,      taumod,     taunorm,     tauref0
    double precision ::   temp(ip11),  toxx(ip12),  toxy(ip12),  toxz(ip12),  toyy(ip12)
    double precision ::   toyz(ip12),  tozz(ip12),           u,       utaur,       utsnu
    double precision ::          utx,        utxt,         uty,        utyt,         utz
    double precision ::         utzt,           v,           w,     x(ip21),       xcfac
    double precision ::      y(ip21),       ycfac,     z(ip21),       zcfac
    logical          :: ecrcom
!
!-----------------------------------------------------------------------
!
    character(len=1 ) :: c
!
!
!     double cote
    c=char(34)
!
    kimpl=kimp
    kimpl=1
!
!     -------------------------------------------------------
!     ecriture commemtaires dans "fsor2"
    ecrcom=.true.
!
!     -------------------------------------------------------
!     adimensionnement du frottement
!
!     definitions dans "dfph.f" :
!     gam=GAMMA EAU ; gam1=gam-1 ; gam2=0.5*gam1
!
!     Coefficient de frottement parietal ramene aux grandeurs a l'infini amont:
!     Cf0=tau_p / 0.5 rho_inf U_inf**2  : grandeurs avec dimensions
!
!     Le tenseur des contraintes est adimensionne par:
!     rho_i  -> masse volumique conditions d'arret
!     a_i**2 -> vitesse du son dans les conditions d'arret au carre
!
!     Cf0   =   (tau_calcul) * tauref
!
!     tauref = 2 * rho_i * a_i**2 / ( rho_inf * U_inf**2 )
!                  <-----grandeurs avec dimension------>
!
    tauref0=2.*(1.+gam2*rm0**2)**(gam/gam1) / rm0**2
!
!---------------------------------------------------------------
!     adimensionnement du flux de chaleur parietal
!
!     coefficient de chaleur parietal ramene aux grandeurs a l'infini amont:
!     Ch= phi / (rho_inf U_inf hi_inf) : grandeurs avec dimension
!
!     le flux de chaleur parietal phi est adimensionne par: rho_i *a_i**3
!
!     Ch  =  (phi_calcul) * phiref0
!
!     phiref0 = (rho_i a_i**3)/(rho_inf U_inf hi_inf) =
!             = (gamma-1)*(rho_i a_i)/(rho_inf a_inf M_inf)
!              <-----grandeurs avec dimension-------->
!
    phiref0=( gam1*(1.+gam2*rm0**2)**(1./gam1+0.5) )/rm0
!
!     -------------------------------------------------------
!     SORTIES RELATIVES A DES VALEURS SUR LES PAROIS
!
    if(rank==0)then
    open(sor2 ,file='pres')
    write(sor2,'(''TITLE='',a1,a80,a1)')c,titrt1,c
    if(ecrcom) then
!        write(sor2,'("#==>utsorfr: adimensionnement du frottement ",
!     &  "cf0=0.5 rho_inf U_inf**2",/,
!     &  "#  gamma=",1pe11.3,5x,"rm0=",1pe11.3,2x,"dans fdon1",4x,
!     &  "tau_ref_0=",1pe11.3,3x,"phiref=",1pe11.3)')gam,rm0,tauref0,
!     &  phiref0
!        write(sor2,'("#",/,"#",t5,"x/L",t15,"y/L",t25,"z/L",
!     &  t33,"t_x/RU2_0",t44,"t_y/RU2_0",t55,"t_z/RU2_0",t67,"Cf0",
!     &  t78,"Kp",t89,"p/pi",t99,"Mp_isent",t114,"qx",t125,"qy",t136,
!     &  "qz",t144,"q_tot",t155,"m1   m2",t164,"L Utau/nu  U_tau/a_i")')
!
       write(sor2,'(''VARIABLES = '',a1,17(a,a1,'', '',a1),a,a1)') &
            c,'x/L',c, c,'y/L',c, c,'z/L',c, c,'t_x/rU2',c, c,'t_y/rU2',c, &
            c,'t_z/rU2',c, c,'Cf0',c, c,'Kp',c, c,'P/Pi',c, &
            c,'Mp_isen',c, c,'qx',c, c,'qy',c, c,'qz',c, &
            c,'i',c, c,'j',c, c,'L Utau/nu',c, c,'Utau/a_i',c
       write(sor2,'("# Mach0=",1pe11.3,4x,"tau_ref_0=",1pe11.3,' &
            //'3x,"phiref=",1pe11.3)')rm0,tauref0,phiref0
       write(sor2,'("#",t5,"x/L",t15,"y/L",t25,"z/L",' &
            //'t33,"t_x/rU2",t44,"t_y/rU2",t55,"t_z/rU2",t67,"Cf0",' &
            //'t78,"Kp",t89,"P/Pi",t99,"Mp_isen",t114,"qx",t125,"qy",t136,' &
            //'"qz",t144,"i   j",t164,"L Utau/nu",t175,' &
            //'"Utau/a_i")')
    close(sor2)
    end if
    end if
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
    if(abs(v0).le.tiny(1.)) then
       write(imp,'(/,"!!!utsorfr: v0=0  devient 1")')
       v0=1.
    end if
    if(abs(q0spi0).le.tiny(1.)) then
       write(imp,'(/,"!!!utsorfr: q0spi0=0  devient 1")')
       q0spi0=1.
    end if
    if(abs(pa1).le.tiny(1.)) then
       write(imp,'(/,"!!!utsorfr: pa1=0  devient 1")')
       pa1=1.
    end if
!
!      write(imp,'("===>utsorfr: nbfll=",i3)')nbfll
!
!    if(nbfll.eq.0) then !pas de paroi a traiter
!       return
!    endif
!
    do mf=1,nbfll
!       boucle sur les parois
!
       mfl=nmfint(mf)
       mfg=bcint_to_bcintg(mf)
       call start_keep_order(mfg,bcintg_to_proc)
       open(sor2 ,file='pres',position="append")
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
       if(kimpl.eq.1) then
          write(imp,987) bl_to_bg(l),bcl_to_bcg(mfl),iminf,imaxf,jminf,jmaxf,kminf,kmaxf
987       format('integration des pressions par frontiere :', &
               /1x,39('-') &
               //1x,'zone ',i3,'  - frontiere ',i3/1x,26('-'), &
               //5x,'imin = ',i3,5x,'imax = ',i3, &
               //5x,'jmin = ',i3,5x,'jmax = ',i3, &
               //5x,'kmin = ',i3,5x,'kmax = ',i3/)
       endif
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
!       initialisation pression
       cxav=0.
       cyav=0.
       czav=0.
       clav=0.
       cmav=0.
       cnav=0.
       syz =0.
       sxz =0.
       sxy =0.
       sml =0.
!       initialisation frottement
       cxavfr=0.
       cyavfr=0.
       czavfr=0.
       clavfr=0.
       cmavfr=0.
       cnavfr=0.
       syzfr =0.
       sxzfr =0.
       sxyfr =0.
       smlfr =0.
!
       do m2=m2min,m2maxm1
!         boucle sur les bandes
!
!         initialisation pression
          cxavb=0.
          cyavb=0.
          czavb=0.
          clavb=0.
          cmavb=0.
          cnavb=0.
!
!         initialisation frottement
          cxavbfr=0.
          cyavbfr=0.
          czavbfr=0.
          clavbfr=0.
          cmavbfr=0.
          cnavbfr=0.
!
!         surface de la bande bande
          sxyb=0.
!
!                 mfac   : pointeur tableaux toutes frontieres
!                 mfacn  : pointeur tableaux frontieres a normales stockees
!                 nfacf  : pointeur tableaux toutes cellules
!
          do m1=m1min,m1maxm1
!           boucle sur les cellules de la bande
             mfac =m0b+m1+(m2-1)*idm
             mfacn=m0n+m1+(m2-1)*idm
             nfacf=ncbd(mfac)
             nfac1=nfacf-n0c+n0n+idfac
             nfac2=nfacf-n0c+n0n+idfac+idf1
             nfac3=nfacf-n0c+n0n+idfac+idf1+idf2
             nfac4=nfacf-n0c+n0n+idfac+idf2
             nn   =nfacf-n0c+n0n
!
!           pression statique
             pspi0=ps(nfacf)/pa1
             xcfac=(x(nfac1)+x(nfac2)+x(nfac3)+x(nfac4))/4.
             ycfac=(y(nfac1)+y(nfac2)+y(nfac3)+y(nfac4))/4.
             zcfac=(z(nfac1)+z(nfac2)+z(nfac3)+z(nfac4))/4.
             dx1  =x(nfac3)-x(nfac1)
             dy1  =y(nfac3)-y(nfac1)
             dz1  =z(nfac3)-z(nfac1)
             dx2  =x(nfac4)-x(nfac2)
             dy2  =y(nfac4)-y(nfac2)
             dz2  =z(nfac4)-z(nfac2)
             dsyz =abs(dy1*dz2-dz1*dy2)/2.
             dsxz =abs(dz1*dx2-dx1*dz2)/2.
             dsxy =abs(dx1*dy2-dy1*dx2)/2.
             dsml =sqrt(dsyz*dsyz+dsxz*dsxz+dsxy*dsxy)
             dcx  =(p0spi0-pspi0)*dsyz*sign(1.,nxn(mfacn))
             dcy  =(p0spi0-pspi0)*dsxz*sign(1.,nyn(mfacn))
             dcz  =(p0spi0-pspi0)*dsxy*sign(1.,nzn(mfacn))
             dcl  =dcy*(zcfac-zref)-dcz*(ycfac-yref)
             dcm  =dcx*(zcfac-zref)-dcz*(xcfac-xref)
             dcn  =dcx*(ycfac-yref)-dcy*(xcfac-xref)
             cxav =cxav+dcx
             cyav =cyav+dcy
             czav =czav+dcz
             clav =clav+dcl
             cmav =cmav+dcm
             cnav =cnav+dcn
             syz  =syz+dsyz
             sxz  =sxz+dsxz
             sxy  =sxy+dsxy
             sml  =sml+dsml
             cxavb=cxavb+dcx
             cyavb=cyavb+dcy
             czavb=czavb+dcz
             clavb=clavb+dcl
             cmavb=cmavb+dcm
             cnavb=cnavb+dcn
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
             utxt   =utx-taunorm*nxn(mfacn)
             utyt   =uty-taunorm*nyn(mfacn)
             utzt   =utz-taunorm*nzn(mfacn)
             dcxfr  =utxt*dsml
             dcyfr  =utyt*dsml
             dczfr  =utzt*dsml
             dclfr  =dcyfr*(zcfac-zref)-dczfr*(ycfac-yref)
             dcmfr  =dcxfr*(zcfac-zref)-dczfr*(xcfac-xref)
             dcnfr  =dcxfr*(ycfac-yref)-dcyfr*(xcfac-xref)
             cxavfr =cxavfr+dcxfr
             cyavfr =cyavfr+dcyfr
             czavfr =czavfr+dczfr
             clavfr =clavfr+dclfr
             cmavfr =cmavfr+dcmfr
             cnavfr =cnavfr+dcnfr
             cxavbfr=cxavbfr+dcxfr
             cyavbfr=cyavbfr+dcyfr
             czavbfr=czavbfr+dczfr
             clavbfr=clavbfr+dclfr
             cmavbfr=cmavbfr+dcmfr
             cnavbfr=cnavbfr+dcnfr
!
!           surface
             sxyb=sxyb+dsxy
!           fin de boucle sur les cellules de la bande
          enddo
!
          if(sxyb.gt.0.) then
!
!           bande de surface non nulle
!           pression
             cxavb  =cxavb/(q0spi0*sxyb)
             cyavb  =cyavb/(q0spi0*sxyb)
             czavb  =czavb/(q0spi0*sxyb)
             clavb  =clavb/(q0spi0*sxyb*xlref)
             cmavb  =cmavb/(q0spi0*sxyb*xlref)
             cnavb  =cnavb/(q0spi0*sxyb*xlref)
             cxaerob= cxavb*csal*csbe-cyavb*snbe+czavb*snal*csbe
             cyaerob= cxavb*csal*snbe+cyavb*csbe+czavb*snal*snbe
             czaerob=-cxavb*snal+czavb*csal
             claerob= clavb*csal*csbe+cmavb*snbe+cnavb*snal*csbe
             cmaerob=-clavb*csal*snbe+cmavb*csbe-cnavb*snal*snbe
             cnaerob=-clavb*snal+cnavb*csal
!
!           frottement
             cxavbfr=cxavbfr*tauref0/ sxyb
             cyavbfr=cyavbfr*tauref0/ sxyb
             czavbfr=czavbfr*tauref0/ sxyb
             clavbfr=clavbfr*tauref0/(sxyb*xlref)
             cmavbfr=cmavbfr*tauref0/(sxyb*xlref)
             cnavbfr=cnavbfr*tauref0/(sxyb*xlref)
             cxaerfr= cxavbfr*csal*csbe-cyavbfr*snbe+czavbfr*snal*csbe
             cyaerfr= cxavbfr*csal*snbe+cyavbfr*csbe+czavbfr*snal*snbe
             czaerfr=-cxavbfr*snal+czavbfr*csal
             claerfr= clavbfr*csal*csbe+cmavbfr*snbe+cnavbfr*snal*csbe
             cmaerfr=-clavbfr*csal*snbe+cmavbfr*csbe-cnavbfr*snal*snbe
             cnaerfr=-clavbfr*snal+cnavbfr*csal
!
             if(kimpl.eq.1) then
                write(imp,991) m2,sxyb,cxavb,cyavb,czavb,clavb,cmavb, &
                     cnavb,cxaerob,cyaerob,czaerob,claerob,cmaerob,cnaerob
991             format(//,1x,'bande numero ',i3,' :',/,1x,('-'),// &
                     /,1x,'surface mouillee projetee sur xy : ',e12.4,/ &
                     /,1x,'efforts pression dans le repere avion : ',/ &
                     /,5x,'cx = ',f8.4,5x,'cy = ',f8.4,5x,'cz = ',f8.4 &
                     ,5x,'cl = ',f8.4,5x,'cm = ',f8.4,5x,'cn = ',f8.4,// &
                     /,1x,'efforts pression dans le repere aerodynamique : ',/ &
                     /,5x,'cx = ',f8.4,5x,'cy = ',f8.4,5x,'cz = ',f8.4 &
                     ,5x,'cl = ',f8.4,5x,'cm = ',f8.4,5x,'cn = ',f8.4//)
!
                write(imp,881) cxavbfr,cyavbfr,czavbfr, &
                     clavbfr,cmavbfr,cnavbfr, cxaerfr,cyaerfr,czaerfr, &
                     claerfr,cmaerfr,cnaerfr
881             format( &
                     /,1x,'efforts frottement dans le repere avion : ',/ &
                     /,5x,'cx = ',f8.4,5x,'cy = ',f8.4,5x,'cz = ',f8.4 &
                     ,5x,'cl = ',f8.4,5x,'cm = ',f8.4,5x,'cn = ',f8.4,// &
                     /,1x,'efforts pression dans le repere aerodynamique : ',/ &
                     /,5x,'cx = ',f8.4,5x,'cy = ',f8.4,5x,'cz = ',f8.4 &
                     ,5x,'cl = ',f8.4,5x,'cm = ',f8.4,5x,'cn = ',f8.4//)
             endif
!
          endif
!
!         ecriture des variables aux centres des mailles surfaciques :
!
!         x, y, z, cfx, cfy, cfz, |cf0|, kp, ps
!
          do m1=m1min,m1maxm1
!           boucle sur les facettes de la bande
             mfac =m0b+m1+(m2-1)*idm
             nfacm=m0n+m1+(m2-1)*idm
             nfacf=ncbd(mfac)
             nfac1=nfacf-n0c+n0n+idfac
             nfac2=nfacf-n0c+n0n+idfac+idf1
             nfac3=nfacf-n0c+n0n+idfac+idf1+idf2
             nfac4=nfacf-n0c+n0n+idfac+idf2
             nn   =nfacf-n0c+n0n
!
!           pression statique / pi(infini amont)
             pspi0 =ps(nfacf)/pa1
             xcfac =(x(nfac1)+x(nfac2)+x(nfac3)+x(nfac4))/4.
             ycfac =(y(nfac1)+y(nfac2)+y(nfac3)+y(nfac4))/4.
             zcfac =(z(nfac1)+z(nfac2)+z(nfac3)+z(nfac4))/4.
             u     =s(nfacf,2)/(s(nfacf,1)*v0)
             v     =s(nfacf,3)/(s(nfacf,1)*v0)
             w     =s(nfacf,4)/(s(nfacf,1)*v0)
             qq    =(u*u+v*v+w*w)*v0*v0
             rm =sqrt(qq)/cson(nfacf)
!           pression d arret
!            pa    =ps*(1.+gam2*rm2)**(gam*gam4)
!           pression arret / pi(infini amont)
!            paspa1=pa/pa1
             akp=(pspi0-p0spi0)/q0spi0
!
!           Mach parietal isentropique
!           p/pi = ( 1+ (gamma-1)/2 M^2 ) ^ -gamma/(gamma-1)
!
!            xmeisent2=(pspi0**(-(gam-1.)/gam)-1.)/
!     &                 (0.5*(gam-1.))
!            xmeisent =sqrt(max(0.,xmeisent2))
!
!           frottement
             utx=toxx(nfacf)*nxn(nfacm)+toxy(nfacf)*nyn(nfacm)+ &
                  toxz(nfacf)*nzn(nfacm)
             uty=toxy(nfacf)*nxn(nfacm)+toyy(nfacf)*nyn(nfacm)+ &
                  toyz(nfacf)*nzn(nfacm)
             utz=toxz(nfacf)*nxn(nfacm)+toyz(nfacf)*nyn(nfacm)+ &
                  tozz(nfacf)*nzn(nfacm)
!           projection du frottement sur la surface
             taunorm=utx*nxn(nfacm)+uty*nyn(nfacm)+utz*nzn(nfacm)
             utxt   =(utx-taunorm*nxn(nfacm))*tauref0
             utyt   =(uty-taunorm*nyn(nfacm))*tauref0
             utzt   =(utz-taunorm*nzn(nfacm))*tauref0
             cfinf0 =sqrt(utxt**2+utyt**2+utzt**2)
             cfinf0 =cfinf0*sign(1.,utxt)
!
!           flux de chaleur
             phipadi=qcx(nfacf)*nxn(nfacm)+qcy(nfacf)*nyn(nfacm)+ &
                  qcz(nfacf)*nzn(nfacm)
!EG d       composante tangentielle de q
             qtx=qcx(nfacf)-phipadi*nxn(nfacm)
             qty=qcy(nfacf)-phipadi*nyn(nfacm)
             qtz=qcz(nfacf)-phipadi*nzn(nfacm)
             qt=max(1.e-10,sqrt(qtx**2+qty**2+qtx**2))
             t1=qtx/qt
             t2=qty/qt
             t3=qtz/qt
             phit=(qcx(nfacf)*t1+qcy(nfacf)*t2+qcz(nfacf)*t3)
             phip=phipadi
!            phip   =phipadi*phiref0
!            phip=sqrt(qcx(nfacf)**2+qcy(nfacf)**2+qcz(nfacf)**2)*
!     &               phiref0
!EG f
!
!           utaur : U_tau / a_i
!           utsnu : Lref * u_tau / nu
             nfaci=ncin(mfac)
             taumod =sqrt(utxt**2+utyt**2+utzt**2)/tauref0
             utaur=sqrt(taumod/s(nfacf,1))
             utsnu=sqrt(s(nfacf,1)*taumod)/mu(nfaci)
!
             write(sor2,555) xcfac,ycfac,zcfac,utxt,utyt,utzt,cfinf0, &
                  akp,ps(nfacf),qcx(nfacf),qcy(nfacf),phit,phip, &
                  m1,m2,utsnu,utaur
555          format(3(1pe12.4),6(1pe11.3),2x,4(1pe11.3),2i4,2(1pe11.3))
!           fin de boucle sur les facettes de la bande
          enddo
!
!         fin de boucle sur les bandes
       enddo

!
!       pression
       cxav   =cxav/(q0spi0*sref)
       cyav   =cyav/(q0spi0*sref)
       czav   =czav/(q0spi0*sref)
       clav   =clav/(q0spi0*sref*xlref)
       cmav   =cmav/(q0spi0*sref*xlref)
       cnav   =cnav/(q0spi0*sref*xlref)
       cxaero = cxav*csal*csbe-cyav*snbe+czav*snal*csbe
       cyaero = cxav*csal*snbe+cyav*csbe+czav*snal*snbe
       czaero =-cxav*snal+czav*csal
       claero = clav*csal*csbe+cmav*snbe+cnav*snal*csbe
       cmaero =-clav*csal*snbe+cmav*csbe-cnav*snal*snbe
       cnaero =-clav*snal+cnav*csal
       cxavtot=cxavtot+cxav
       cyavtot=cyavtot+cyav
       czavtot=czavtot+czav
       clavtot=clavtot+clav
       cmavtot=cmavtot+cmav
       cnavtot=cnavtot+cnav
!
!       frottement
       cxavfr =cxavfr*tauref0/ sref
       cyavfr =cyavfr*tauref0/ sref
       czavfr =czavfr*tauref0/ sref
       clavfr =clavfr*tauref0/(sref*xlref)
       cmavfr =cmavfr*tauref0/(sref*xlref)
       cnavfr =cnavfr*tauref0/(sref*xlref)
       cxaefr = cxavfr*csal*csbe-cyavfr*snbe+czavfr*snal*csbe
       cyaefr = cxavfr*csal*snbe+cyavfr*csbe+czavfr*snal*snbe
       czaefr =-cxavfr*snal+czavfr*csal
       claefr = clavfr*csal*csbe+cmavfr*snbe+cnavfr*snal*csbe
       cmaefr =-clavfr*csal*snbe+cmavfr*csbe-cnavfr*snal*snbe
       cnaefr =-clavfr*snal+cnavfr*csal
       cxavtfr=cxavtfr+cxavfr
       cyavtfr=cyavtfr+cyavfr
       czavtfr=czavtfr+czavfr
       clavtfr=clavtfr+clavfr
       cmavtfr=cmavtfr+cmavfr
       cnavtfr=cnavtfr+cnavfr
!
       if(kimpl.eq.1) then
          write(imp,985)
          write(imp,989) sml,sxy,syz,sxz,cxav,cyav,czav,clav,cmav, &
               cnav,cxaero,cyaero,czaero,claero,cmaero,cnaero
985       format(//,1x,'frontiere complete :',/,1x,20('-'),/)
989       format(/,1x,'surface mouillee          : ',e12.4,/ &
               /,1x,'surface mouillee projetee sur xy : ',e12.4,/ &
               /,1x,'surface mouillee projetee sur yz : ',e12.4,/ &
               /,1x,'surface mouillee projetee sur xz : ',e12.4,// &
               /,1x,'efforts pression dans le repere avion : ',/ &
               /,5x,'cx = ',f8.4,5x,'cy = ',f8.4,5x,'cz = ',f8.4 &
               ,5x,'cl = ',f8.4,5x,'cm = ',f8.4,5x,'cn = ',f8.4,// &
               /,1x,'efforts pression dans le repere aerodynamique : ',/ &
               /,5x,'cx = ',f8.4,5x,'cy = ',f8.4,5x,'cz = ',f8.4 &
               ,5x,'cl = ',f8.4,5x,'cm = ',f8.4,5x,'cn = ',f8.4)
!
          write(imp,889) cxavfr,cyavfr,czavfr,clavfr,cmavfr, &
               cnavfr, cxaefr,cyaefr,czaefr,claefr,cmaefr,cnaefr
889       format( &
               /,1x,'efforts frottement dans le repere avion : ',/ &
               /,5x,'cx = ',f8.4,5x,'cy = ',f8.4,5x,'cz = ',f8.4 &
               ,5x,'cl = ',f8.4,5x,'cm = ',f8.4,5x,'cn = ',f8.4,// &
               /,1x,'efforts frottement dans le repere aerodynamique : ',/ &
               /,5x,'cx = ',f8.4,5x,'cy = ',f8.4,5x,'cz = ',f8.4 &
               ,5x,'cl = ',f8.4,5x,'cm = ',f8.4,5x,'cn = ',f8.4)
       endif
       close(sor2)
      call END_KEEP_ORDER(mfg,bcintg_to_proc)
!       fin de boucle sur les parois
    enddo
!
    call SUM_MPI(cxavtot)
    call SUM_MPI(cyavtot)
    call SUM_MPI(czavtot)
    call SUM_MPI(clavtot)
    call SUM_MPI(cmavtot)
    call SUM_MPI(cnavtot)
!     pression
    cxaero= cxavtot*csal*csbe-cyavtot*snbe+czavtot*snal*csbe
    cyaero= cxavtot*csal*snbe+cyavtot*csbe+czavtot*snal*snbe
    czaero=-cxavtot*snal+czavtot*csal
    claero= clavtot*csal*csbe+cmavtot*snbe+cnavtot*snal*csbe
    cmaero=-clavtot*csal*snbe+cmavtot*csbe-cnavtot*snal*snbe
    cnaero=-clavtot*snal+cnavtot*csal
!
    call SUM_MPI(cxavtfr)
    call SUM_MPI(cyavtfr)
    call SUM_MPI(czavtfr)
    call SUM_MPI(clavtfr)
    call SUM_MPI(cmavtfr)
    call SUM_MPI(cnavtfr)
!     frottement
    cxaetfr= cxavtfr*csal*csbe-cyavtfr*snbe+czavtfr*snal*csbe
    cyaetfr= cxavtfr*csal*snbe+cyavtfr*csbe+czavtfr*snal*snbe
    czaetfr=-cxavtfr*snal+czavtfr*csal
    claetfr= clavtfr*csal*csbe+cmavtfr*snbe+cnavtfr*snal*csbe
    cmaetfr=-clavtfr*csal*snbe+cmavtfr*csbe-cnavtfr*snal*snbe
    cnaetfr=-clavtfr*snal+cnavtfr*csal
!
    if(kimpl.eq.1.and.rank==0) then
       write(imp,890) cxavtot,cyavtot,czavtot,clavtot,cmavtot,cnavtot, &
            cxaero,cyaero,czaero,claero,cmaero,cnaero
890    format('efforts globaux pression - configuration ', &
            'complete :'/,1x,40('-'),//,1x,'repere avion :',/ &
            /,5x,'cx = ',f8.4,5x,'cy = ',f8.4,5x,'cz = ',f8.4 &
            ,5x,'cl = ',f8.4,5x,'cm = ',f8.4,5x,'cn = ',f8.4,// &
            /,1x,'repere aerodynamique :',/ &
            /,5x,'cx = ',f8.4,5x,'cy = ',f8.4,5x,'cz = ',f8.4 &
            ,5x,'cl = ',f8.4,5x,'cm = ',f8.4,5x,'cn = ',f8.4)
!
       write(imp,990) cxavtfr,cyavtfr,czavtfr,clavtfr,cmavtfr,cnavtfr , &
            cxaetfr,cyaetfr,czaetfr,claetfr,cmaetfr,cnaetfr
990    format( &
            /,1x,'efforts globaux frottement repere avion :',/ &
            /,5x,'fx = ',1pe11.3,5x,'fy = ',1pe11.3,5x,'fz = ',1pe11.3 &
            ,5x,'mx = ',1pe11.3,5x,'cm = ',1pe11.3,5x,'cn = ',1pe11.3,// &
            /,1x,'repere aerodynamique :',/ &
            /,5x,'cx = ',f8.4,5x,'cy = ',f8.4,5x,'cz = ',f8.4 &
            ,5x,'cl = ',f8.4,5x,'cm = ',f8.4,5x,'cn = ',f8.4)
    endif
!
    return
  end subroutine utsorfr
end module mod_utsorfr
