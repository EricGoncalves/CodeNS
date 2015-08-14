module mod_met_intep3
  implicit none
contains
  subroutine met_intep3( &
       ncbd,ncin,s, &
       sn,vol, &
       dist,mnpar,mu, &
       pi,tau,us,ut,un, &
       toxx,toxy,toxz,toyy,toyz,tozz, &
       x,y,z,nxn,nyn,nzn, &
       ps,cson,temp, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des epaisseurs integrales de la couche limite.
!_A    Boucle sur les parois et double boucle sur les facettes.
!_A    ecriture des epaisseurs dans le fichier "fateps" etiquette sor3
!      integration jusqu'a delta vrai et non pas jusqu'au noeud superieur
!_A
!
!     INP
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    s          : arg real(ip11,ip60 ) ; variables de calcul
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_L    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_L                                        contraintes visqueuses
!_L    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_L                                        contraintes visqueuses
!_L    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_L                                        contraintes visqueuses
!_L    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_L                                        contraintes visqueuses
!_L    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_L                                        contraintes visqueuses
!_L    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_L                                        contraintes visqueuses
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
!_I    gam        : com real             ; rapport des chaleurs specifiques
!_I    gam1       : com real             ; rap chal spec -1
!_I    gam2       : com real             ; (rap chal spec -1)/2
!_I    gam4       : com real             ; 1/(rap chal spec -1)
!_I    cp         : com real             ; chal spec a pres cste adim
!_I    cv         : com real             ; chal spec a vol cst adim
!_I    pr         : com real             ; nombre de Prandtl
!_I    prt        : com real             ; nombre de Prandtl turbulent
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    kvglo      : com int              ; cle calcul dea grandeurs globales
!_I    nbfll      : com int              ; nb de front d'integration des
!_I                                        grandeurs globales
!_I    nmfint     : com int (mtb       ) ; no des front d'integration des
!_I                                        grandeurs globales
!
!     OUT
!
!     I/O
!
!     LOC
!_L    pi         : arg real(ip00      ) ; pression d'arret locale
!_L    ps         : arg real(ip00      ) ; pression statique locale
!_L    tau        : arg real(ip00      ) ; frottement local total
!_L    us         : arg real(ip00      ) ; vitesse suivant S ou tab de travail
!_L    ut         : arg real(ip00      ) ; vitesse suivant T
!_L    un         : arg real(ip00      ) ; vitesse suivant N dans
!_L                                        repere couche limite
!_I    dvxx       : arg real(ip00      ) ; composante xx du gradient de vitesse
!_I    dvxy       : arg real(ip00      ) ; composante xy du gradient de vitesse
!_I    dvxz       : arg real(ip00      ) ; composante xz du gradient de vitesse
!_I    dvyx       : arg real(ip00      ) ; composante yx du gradient de vitesse
!_I    dvyy       : arg real(ip00      ) ; composante yy du gradient de vitesse
!_I    dvyz       : arg real(ip00      ) ; composante yz du gradient de vitesse
!_I    dvzx       : arg real(ip00      ) ; composante zx du gradient de vitesse
!_I    dvzy       : arg real(ip00      ) ; composante zy du gradient de vitesse
!_I    dvzz       : arg real(ip00      ) ; composante zz du gradient de vitesse
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use proprieteflu
    use chainecarac
    use definition
    use sortiefichier
    use constantes
    use modeleturb
    use mod_teq_gradv
    implicit none
    integer          ::          i1,         i2,       idd1,       idd2,       idd3
    integer          ::        ideb,       idf1,       idf2,      idfac,        idm
    integer          ::        idm3,       ierr,         ii,      imaxf,      iminf
    integer          ::      isens3,         j1,         j2,       jdd1,       jdd2
    integer          ::        jdd3,       jdeb,         jj,      jmaxf,      jminf
    integer          ::          k1,         k2,       kdd1,       kdd2,       kdd3
    integer          ::        kdeb,         kk,      kmaxf,      kminf,          l
    integer          ::        lbgr,          m,        m0b,        m0n,     m0ndeb
    integer          ::      m0nfin,         m1,      m1max,    m1maxm1,      m1min
    integer          ::          m2,      m2max,    m2maxm1,      m2min,         m3
    integer          ::       m3del,     m3delp,     m3delt,     m3delv,      m3max
    integer          ::       m3min,      m3mxd,      m3pmx,      m3tmx,      m3vmx
    integer          ::        mdel,      mdel2,         mf,       mfac,      mfacn
    integer          ::         mfl,mnpar(ip12),       mpar,          n,        n0c
    integer          ::         n0n,        nc0, ncbd(ip41),        nci, ncin(ip41)
    integer          ::         ncj,        nck,       ndel,      ndel2,      ndelp
    integer          ::       ndelt,      ndelv,      nfac1,      nfac2,      nfac3
    integer          ::       nfac4,      nfacf,        nid,       nijd,        njd
    double precision ::           am2,  cmui1(ip21),  cmui2(ip21),  cmuj1(ip21),  cmuj2(ip21)
    double precision ::   cmuk1(ip21),  cmuk2(ip21),          csc,   cson(ip11),        ddist
    double precision ::          del1,        del1i,          dev,   dist(ip12),       distm1
    double precision ::             e,        epspi,       epstau,       epsvor,         hpar
    double precision ::         hpari,     mu(ip12),           nx,    nxn(ip42),           ny
    double precision ::     nyn(ip42),           nz,    nzn(ip42),            p,     pi(ip00)
    double precision ::         pimax,     ps(ip11),           qq,         reyl,         rhoe
    double precision ::           rm2,       rm3dlp,       rm3dlt,       rm3dlv,           ro
    double precision ::         rpdel,           rr,         rus0,         rus1, s(ip11,ip60)
    double precision ::            sg,sn(ip31*ndir),        snorm,         somd,        somru
    double precision ::        somru2,         somu,        somu2,           sx,           sy
    double precision ::            sz,    tau(ip00),        taumx,           te,   temp(ip11)
    double precision ::       theta11,     theta11i,         tmod,   toxx(ip12),   toxy(ip12)
    double precision ::    toxz(ip12),   toyy(ip12),   toyz(ip12),   tozz(ip12),          tue
    double precision ::            tx,           ty,           tz,            u,        uemod
    double precision ::           uen,          uex,         uex1,         uex2,          uey
    double precision ::          uey1,         uey2,          uez,         uez1,         uez2
    double precision ::      un(ip00),          und,     us(ip00),          us0,          us1
    double precision ::           usd,     ut(ip00),          utd,          utx,          uty
    double precision ::           utz,            v,    vol(ip11),        vormx,           vv
    double precision ::             w,      x(ip12),        xcfac,          xme,         xmue
    double precision ::            xn,           xs,           xt,      y(ip12),        ycfac
    double precision ::            yn,           ys,           yt,      z(ip12),        zcfac
    double precision ::            zn,           zs,           zt
    logical          ::    iok,ouvert
    double precision,allocatable :: dvxx(:),dvxy(:),dvxz(:),dvyx(:),dvyy(:)
    double precision,allocatable :: dvyz(:),dvzx(:),dvzy(:),dvzz(:),vort(:)
!
!-----------------------------------------------------------------------
!
    character(len=20) ::  nomfic
    character(len=1 ) :: c
    character(len=5 ) :: control
!
    ALLOCATE(dvxx(ip00),dvxy(ip00),dvxz(ip00),dvyx(ip00),dvyy(ip00),dvyz(ip00), &
         dvzx(ip00),dvzy(ip00),dvzz(ip00),vort(ip12))


!     double cote
    c=char(34)
!
!     -------------------------------------------------------
!
!      write(imp,'(/,"===>met_intep3: nbfll=",i3)')nbfll
!
    if(nbfll.eq.0) then
!       pas de paroi a traiter
       write(imp,'("!!!!met_intep3: pas de paroi a traiter")')
       return
    endif
!
    if(kcaldis.eq.0 .and. klecdis.eq.0) then
!       pas de calcul de distances
       write(imp,'("!!!!met_intep3: pas de calcul de distances")')
       return
    end if
!
    if(epspid.le.0. .or. epstaud.le.0. .or. epsvord.le.0) then
!       pas de definition des constantes pour delta
       write(imp,'("!!!!met_intep: pas de definition des constantes dans fatdon")')
       return
    endif
!
!     --------------------------------------------------------
!     ouverture fichier
!
!     inquire(sor3,opened=ouvert,name=nomfic)
!      if(ouvert) then
!        write(imp,'("!!!!met_intep: fichier fateps etiquette=",i3,
!     &  2x,"deja ouvert. ",/,12x,"nom=",a20,/,12x,"RETURN")')
!     &  sor3,nomfic
!        return
!      end if
!      open(sor3,file='fateps3',form='formatted',err=100)
!
!     &    c,'Hi',c, c,'rhoe',c, c,'ue_s',c, c,'|Ue|',c, c,'Me',c,
!     &    c,'R_L',c, c,'i',c, c,'j',c, c,'k',c
!
!     ---------------------------------------------------------
!
    lbgr=0
    do mf=1,nbfll
!       boucle sur les parois
!
       mfl=nmfint(mf)
       l=ndlb(mfl)
       if(l.ne.lbgr) then
!         le gradient de la vitesse ne correspond pas au domaine en cours
          lbgr=l
!com      teq_gradv --> grad(v) aux points interieurs au domaine
          call teq_gradv( &
               l, &
               sn, &
               vol,s, &
               ut, &
               dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
               cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
       endif
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
!       m0ndeb,m0nfin : valeurs extremes pointeur facettes paroi
!                       dans tableaux a normales stockees
       m0ndeb=m0n+1
       m0nfin=m0n+mmb(mfl)
!
       ierr=0
       m1min=1
       m2min=1
       if (iminf.eq.imaxf) then
!         frontiere I1 ou I2
!
          m1max=jmaxf-jminf+1
          m1maxm1=m1max-1
          m2max=kmaxf-kminf+1
          m2maxm1=m2max-1
          idf1=ncj
          idf2=nck
!
          if(iminf.eq.i1) then
!           frontiere I1
!
             idfac=nci
             m3min=i1
             m3max=i2-1
             isens3=1
             idm3=nci*isens3
!
          else if(iminf.eq.i2) then
!           frontiere I2
!
             idfac=0
             m3min=i2-1
             m3max=i1
             isens3=-1
             idm3=nci*isens3
          else
             ierr=1
          end if
!         incrementation de "i,j,k" suivant "m1,m2,m3"
!                 ii=ideb+idd1*(m1-1)+idd2*(m2-1)-idd3
!                 jj=jdeb+jdd1*(m1-1)+jdd2*(m2-1)-jdd3
!                 kk=kdeb+kdd1*(m1-1)+kdd2*(m2-1)-kdd3
          idd1=0
          idd2=0
          idd3=isens3
!         jdd1=ncj
          jdd1=1
          jdd2=0
          jdd3=0
          kdd1=0
!         kdd2=nck
          kdd2=1
          kdd3=0
          ideb=m3min
          jdeb=jminf
          kdeb=kminf
!
       elseif (jminf.eq.jmaxf) then
!         frontiere J1 ou J2
!
          m1max=imaxf-iminf+1
          m1maxm1=m1max-1
          m2max=kmaxf-kminf+1
          m2maxm1=m2max-1
          idf1=nci
          idf2=nck
!
          if(jminf.eq.j1) then
!           frontiere J1
!
             idfac=ncj
             m3min=j1
             m3max=j2-1
             isens3=1
             idm3=ncj*isens3
          else if(jminf.eq.j2) then
!           frontiere J2
!
             idfac=0
             m3min=j2-1
             m3max=j1
             isens3=-1
             idm3=ncj*isens3
          else
             ierr=2
          end if
!         incrementation de "i,j,k" suivant "m1,m2,m3"
!                 ii=ideb+idd1*(m1-1)+idd2*(m2-1)-idd3
!                 jj=jdeb+jdd1*(m1-1)+jdd2*(m2-1)-jdd3
!                 kk=kdeb+kdd1*(m1-1)+kdd2*(m2-1)-kdd3
!         idd1=nci
          idd1=1
          idd2=0
          idd3=0
          jdd1=0
          jdd2=0
          jdd3=isens3
          kdd1=0
!         kdd2=nck
          kdd2=1
          kdd3=0
          ideb=iminf
          jdeb=m3min
          kdeb=kminf
!
       elseif (kminf.eq.kmaxf) then
!         frontiere K1 ou K2
!
          m1max=imaxf-iminf+1
          m1maxm1=m1max-1
          m2max=jmaxf-jminf+1
          m2maxm1=m2max-1
          idf1=nci
          idf2=ncj
!
          if(kminf.eq.k1) then
!           frontiere K1
!
             idfac=nck
             m3min=k1
             m3max=k2-1
             isens3=1
             idm3=nck*isens3
          else if(kminf.eq.k2) then
!           frontiere K2
!
             idfac=0
             m3min=k2-1
             m3max=k1
             isens3=-1
             idm3=nck*isens3
!
          else
             ierr=3
          end if
!         incrementation de "i,j,k" suivant "m1,m2,m3"
!                 ii=ideb+idd1*(m1-1)+idd2*(m2-1)-idd3
!                 jj=jdeb+jdd1*(m1-1)+jdd2*(m2-1)-jdd3
!                 kk=kdeb+kdd1*(m1-1)+kdd2*(m2-1)-kdd3
!         idd1=nci
          idd1=1
          idd2=0
          idd3=0
          jdd1=0
!         jdd2=ncj
          jdd2=1
          jdd3=0
          kdd1=0
          kdd2=0
          kdd3=isens3
          ideb=iminf
          jdeb=jminf
          kdeb=m3min
!
       endif
       if(ierr.eq.1) then
          write(imp,'("!!!met_intep: erreur frontiere I1 ou I2")')
          return
       else if(ierr.eq.2) then
          write(imp,'("!!!met_intep: erreur frontiere J1 ou J2")')
          return
       else if(ierr.eq.3) then
          write(imp,'("!!!met_intep: erreur frontiere K1 ou K2")')
          return
       end if
!
       write(sor3,'("ZONE F=POINT, I=",i4," J=",i4)') &
            m1max-m1min,  m2max-m2min

       write(sor3,'("#",t4,"xcfac",t16,"ycfac",t28,"zcfac",' &
            //'t40,"theta11",t54,"H",t64,"theta11i",t77,"Hi",' &
            //'t88,"ue_s",t100,"|Ue|",t111,"Me",t124,"R_L",' &
            //'t135,"i   j   k")')
!
!       ---------------------------------------
!       mfac   : pointeur tableaux toutes frontieres
!       mfacn  : pointeur tableaux frontieres a normales stockees
!       nfacf  : pointeur cellule fictive tableaux toutes cellules
!       nc0    : pointeur cellule adjacente a la paroi tableaux toutes cellules
!       ---------------------------------------
!
       idm=m1max-m1min
       do m2=m2min,m2maxm1
!         boucle sur les bandes
          do m1=m1min,m1maxm1
!           boucle sur les cellules de la bande
!
             mfac =m0b+m1+(m2-1)*idm
             mfacn=m0n+m1+(m2-1)*idm
             nfac1=ncbd(mfac)-n0c+n0n+idfac
             nfac2=ncbd(mfac)-n0c+n0n+idfac+idf1
             nfac3=ncbd(mfac)-n0c+n0n+idfac+idf1+idf2
             nfac4=ncbd(mfac)-n0c+n0n+idfac+idf2
             nfacf=ncbd(mfac)
             nc0  =ncin(mfac)
!           centres des facettes des parois
             xcfac=(x(nfac1)+x(nfac2)+x(nfac3)+x(nfac4))/4.
             ycfac=(y(nfac1)+y(nfac2)+y(nfac3)+y(nfac4))/4.
             zcfac=(z(nfac1)+z(nfac2)+z(nfac3)+z(nfac4))/4.
!           indices "i,j,k" de la cellule adjacente
             ii=ideb+idd1*(m1-1)+idd2*(m2-1)-idd3
             jj=jdeb+jdd1*(m1-1)+jdd2*(m2-1)-jdd3
             kk=kdeb+kdd1*(m1-1)+kdd2*(m2-1)-kdd3
!
             n=nc0-idm3
             do m3=m3min,m3max,isens3
!             boucle suivant "normale" a la paroi-calcul des maximum
!
                ii=ii+idd3
                jj=jj+jdd3
                kk=kk+kdd3
!
                n=n+idm3
                m=n-n0c
                ro=s(n,1)
                u=s(n,2)/ro
                v=s(n,3)/ro
                w=s(n,4)/ro
                e=s(n,5)/ro
                qq=u*u+v*v+w*w
                rhoe=ro*e-.5*ro*qq
!               p=gam1*(rhoe-pinfl)
                p=ps(n)
!               sg=sign(1.,qq)
                rm2=abs(qq*ro/(gam*p))
                rr=y(n)**2+z(n)**2
                vv=qq+omg*(omg*rr+2.*(y(n)*w-z(n)*v))
                sg=sign(1.,vv)
!               am2=abs(vv*ro/(gam*p))
                am2=abs(vv)/cson(n)**2
!
!              Pression d'arret, frottement, vorticite
                pi(m)=sg*abs(p*(1.+gam2*am2)**(gam*gam4)/pa1)
                tau(m)=sqrt(toxx(n)**2+toxy(n)**2+toxz(n)**2+toyy(n)**2+ &
                     toyz(n)**2+tozz(n)**2)
                vort(n)=sqrt( (dvzy(m)-dvyz(m))**2 &
                     +(dvxz(m)-dvzx(m))**2 &
                     +(dvyx(m)-dvxy(m))**2 )
!
             enddo !fin boucle suivant "normale" a la paroi
!
!           --------------------------------------------
!           domaine rattache a la paroi
!           m3mxd : indice "m3" au dela duquel les cellules
!                   de la normale ne sont plus rattachees
!                   a la paroi en cours de calcul
             m3mxd=m3max
             n=nc0-idm3
!
             m3=m3min
             do m3=m3min,m3max,isens3
                n=n+idm3
                if((mnpar(n)-m0ndeb)*(mnpar(n)-m0nfin).gt.0) then
!               cellule non rattachee a la frontiere en cours
                   m3mxd=m3-isens3
                   exit
                endif
             enddo
!
             pimax=-1.
             taumx=-1.
             vormx=-1.
             m3pmx=0
             m3tmx=0
             m3vmx=0
             if(m3mxd.le.0) then
                write(imp,'("!!!met_intep: calcul de maximum impossible m1=",i4,4x,"m2=",i4)')m1,m2
                if(m1.le.2 .and. m2.le.2) then
                   m3=m3min
                   n=nc0-idm3
                   do m3=m3min,m3max,isens3
                      n=n+idm3
                   end do
                endif
                exit
             endif
!
             n=nc0-idm3
             do m3=m3min,m3mxd,isens3
!             boucle suivant "normale" a la paroi,  calcul des maximum
                n=n+idm3
                m=n-n0c
                if(pi(m).gt.pimax) then
                   pimax=pi(m)
                   m3pmx=m3
                end if
                if(tau(m).gt.taumx) then
                   taumx=tau(m)
                   m3tmx=m3
                end if
                if(vort(n).gt.vormx) then
                   vormx=vort(n)
                   m3vmx=m3
                end if
             end do      !fin boucle suivant "normale" a la paroi
!
             if(m3pmx.eq.0) pimax=1.
             if(m3tmx.eq.0) taumx=1.
             if(m3vmx.eq.0) vormx=1.
             if(abs(pimax).le.tiny(1.)) pimax=1.
             if(abs(taumx).le.tiny(1.)) taumx=1.
             if(abs(vormx).le.tiny(1.)) vormx=1.
!
!         --------------------------------------------------------
!           epaisseurs de couche limite.
!           recherche des valeurs seuil pour la vorticite, le frottement
!           total et la pression d'arret au dela du point de frottement
!           maximum, tout en restant dans la zone de rattachement de paroi
!         --------------------------------------------------------
!
             epsvor=epsvord*vormx
             epspi =epspid *pimax
             epstau=epstaud*taumx
             m3delv=0
             m3delp=0
             m3delt=0
             ndelt=ip12
             ndelv=ip12
             ndelp=ip12
!
!           n=nc0-idm3
             n=nc0+(m3tmx-m3min)*isens3*idm3-idm3
             do m3=m3tmx,m3mxd,isens3
!             boucle suivant "normale" a la paroi pour calcul delta
!             recherche de la frontiere au dela du maximum de frottement
!             en n'utilisant que le frottement et le rotationnel
!             delta est compris entre m3delX et m32=m3del+isens3.
!             La cellule m3delX est dans la couche limite
!
                n=n+idm3
                m=n-n0c
!
                if(tau(m).le.epstau .and. m3delt.eq.0) then
                   m3delt=m3-isens3
                   mdel2 =m
                   mdel  =m-idm3
                   ndel2 =n
                   ndelt =n-idm3
                   rm3dlt=(epstau-tau(mdel))/(tau(mdel2)-tau(mdel))
                end if
                if(vort(n).le.epsvor .and. m3delv.eq.0) then
                   m3delv=m3-isens3
                   mdel2 =m
                   mdel  =m-idm3
                   ndel2 =n
                   ndelv =n-idm3
                   if(abs(vort(ndel2)-vort(ndelv)).le.tiny(1.)) then
                      rm3dlv=0.
                   else
                      rm3dlv=(epsvor-vort(ndelv))/(vort(ndel2)-vort(ndelv))
                   endif
                endif
             enddo  !fin boucle suivant "normale" a la paroi pour calcul
!
             n=nc0-idm3
             do m3=m3min,m3mxd,isens3
!             boucle suivant "normale" a la paroi pour calcul delta
!             recherche de la frontiere a partir de la paroi sur Pi uniquement
!
                n=n+idm3
                m=n-n0c
                if(pi(m).ge.epspi .and. m3delp.eq.0) then
                   m3delp=m3-isens3
                   mdel2 =m
                   mdel  =m-idm3
                   ndelp =n-idm3
                   if(abs(pi(mdel2)-pi(mdel)).le.tiny(1.)) then
                      rm3dlp=rm3dlv
                   else
                      rm3dlp=(epspi-pi(mdel))/(pi(mdel2)-pi(mdel))
                   endif
!
                endif
             enddo  !fin boucle suivant "normale" a la paroi pour calcul
!
!           -------------------------------------------------------
!           repere orthonorme de couche limite :
!
!           S  : tangent a la paroi et dans le plan de la vitesse en delta
!           N  : normale a la paroi
!           T  : N vectoriel S. Le repere (S,T,N) est direct
!
!           ->       ->      ->          ->   ->
!           Ue = ues S + uen N  ;  uen = Ue . N
!
!           m3del : indice "m3" de delta
!           ndel  : pointeur cellule de delta (tous domaines)
!                   epaisseur de couche limite prise comme la plus petites
!                   de celles donnees par les differents criteres
!           -------------------------------------------------------
!
             ndel=0
             if(isens3.gt.0) then
                m3del=m3min
                if(m3delt.ne.0) then
                   m3del=m3delt
                   ndel =ndelt
                   rpdel=rm3dlt
                endif
                if(m3delp.ne.0 .and. m3delp.lt.m3del) then
                   m3del=m3delp
                   ndel =ndelp
                   rpdel=rm3dlp
                endif
                if(m3delv.ne.0 .and. m3delv.lt.m3del) then
                   m3del=m3delv
                   ndel =ndelv
                   rpdel=rm3dlv
                endif
             else
                m3del=m3min
                if(m3delt.ne.0) then
                   m3del=m3delt
                   ndel =ndelt
                   rpdel=rm3dlt
                endif
                if(m3delp.ne.0 .and. m3delp.gt.m3del) then
                   m3del=m3delp
                   ndel =ndelp
                   rpdel=rm3dlp
                endif
                if(m3delv.ne.0 .and. m3delv.gt.m3del) then
                   m3del=m3delv
                   ndel =ndelv
                   rpdel=rm3dlv
                endif
             endif
             iok=.false.
!            if(ndel.ne.0) iok=.true.
!            mpar =mnpar(ndel)
!            iok=.false.
             if(ndel.ne.0) mpar =mnpar(ndel)
             if(isens3.gt.0 .and. m3del.gt.m3min) iok=.true.
             if(isens3.lt.0 .and. m3del.lt.m3min) iok=.true.
!
             if(iok) then
!             Calcul de epaisseurs possible. Changement de repere avec
!             les conditions en delta
                uex1=s(ndel,2)/s(ndel,1)
                uey1=s(ndel,3)/s(ndel,1)
                uez1=s(ndel,4)/s(ndel,1)
                ndel2=ndel+idm3
                uex2=s(ndel2,2)/s(ndel2,1)
                uey2=s(ndel2,3)/s(ndel2,1)
                uez2=s(ndel2,4)/s(ndel2,1)
!
                uex=uex1+rpdel*(uex2-uex1)
                uey=uey1+rpdel*(uey2-uey1)
                uez=uez1+rpdel*(uez2-uez1)
                rhoe=s(ndel,1)+rpdel*(s(ndel2,1)-s(ndel,1))
!
                nx =nxn(mpar)
                ny =nyn(mpar)
                nz =nzn(mpar)
                uen=uex*nx+uey*ny+uez*nz
                sx =uex-uen*nx
                sy =uey-uen*ny
                sz =uez-uen*nz
                snorm=sqrt(sx*sx+sy*sy+sz*sz)
                sx=sx/snorm
                sy=sy/snorm
                sz=sz/snorm
!
                tx=ny*sz-nz*sy
                ty=nz*sx-nx*sz
                tz=nx*sy-ny*sx
!
!             inversion de la matrice [S,T,N]
!
                xs=  ty*nz-tz*ny
                xt=-(sy*nz-sz*ny)
                xn=  sy*tz-sz*ty
                ys=-(tx*nz-tz*nx)
                yt=  sy*nz-sz*nx
                yn=-(sx*tz-sz*tx)
                zs=  tx*ny-ty*nx
                zt=-(sx*ny-sy*nx)
                zn=  sx*ty-sy*tx
!
!             passage dans le repere couche limite
!
                n=nc0-idm3
                do m3=m3min,m3max,isens3
!               boucle suivant "normale" a la paroi pour repere couche limite
                   n=n+idm3
                   m=n-n0c
                   us(m)=xs*s(n,2)+ys*s(n,3)+zs*s(n,4)
                   ut(m)=xt*s(n,2)+yt*s(n,3)+zt*s(n,4)
                   un(m)=xn*s(n,2)+yn*s(n,3)+zn*s(n,4)
                enddo  !fin boucle suivant "normale" paroi pour repere cou
!
!             vitesse en delta dans repere couche limite
!
                usd=(xs*uex+ys*uey+zs*uez)*rhoe
                utd=(xt*uex+yt*uey+zt*uez)*rhoe
                und=(xn*uex+yn*uey+zn*uez)*rhoe
!
!             calcul de la deviation beta0
!             frottement - vecteur tension a la paroi T
                utx=toxx(nfacf)*nxn(mfacn)+toxy(nfacf)*nyn(mfacn)+ &
                     toxz(nfacf)*nzn(mfacn)
                uty=toxy(nfacf)*nxn(mfacn)+toyy(nfacf)*nyn(mfacn)+ &
                     toyz(nfacf)*nzn(mfacn)
                utz=toxz(nfacf)*nxn(mfacn)+toyz(nfacf)*nyn(mfacn)+ &
                     tozz(nfacf)*nzn(mfacn)
!             module vecteur tension T
                tmod=sqrt(utx**2+uty**2+utz**2)
!             module vecteur vitesse exterieure Ue
                uemod=sqrt(uex**2+uey**2+uez**2)
!             produit scalaire Ue et T
                tue=utx*uex+uty*uey+utz*uez
!             angle beta0 - deviation dev
!              dev=raddeg*acos(tue/(tmod*uemod))
                dev=0.
!
!             integration des epaisseurs et ecriture sur "sor3"
!
                rus1  =0.
                us1   =0.
                somu  =0.
                somru =0.
                somu2 =0.
                somru2=0.
                somd  =0.
                us0   =0.
                rus0  =0.
                distm1=0.
!
                n=nc0-idm3
                do m3=m3min,m3del,isens3
!               boucle suivant "normale" a la paroi pour integration
!
                   n=n+idm3
                   m=n-n0c
!
                   rus1=us(m)
                   us1=rus1/s(n,1)
                   ddist=dist(n)-distm1
!
                   somu  =somu  +(us0 +us1 )*ddist
                   somru =somru +(rus0+rus1)*ddist
                   somu2 =somu2 +(us0**2  +us1**2  )*ddist
                   somru2=somru2+(rus0*us0+rus1*us1)*ddist
                   somd  =somd  +ddist
!
                   us0   =us1
                   rus0  =rus1
                   distm1=dist(n)
                end do ! fin boucle suivant "normale" a la paroi pour inte
!
!             contribution des integrales entre m3del et delta
                rus1  =usd
                us1   =rus1/rhoe
                ddist =rpdel*(dist(n+idm3)-dist(n))
                somu  =somu  +(us0 +us1 )*ddist
                somru =somru +(rus0+rus1)*ddist
                somu2 =somu2 +(us0**2  +us1**2  )*ddist
                somru2=somru2+(rus0*us0+rus1*us1)*ddist
                somd  =somd  +ddist
!
                somu    =somu /(2.*us1)
                somru   =somru/(2.*rus1)
                somu2   =somu2/(2.*us1**2)
                somru2  =somru2/(2.*rus1*us1)
                del1    =somd-somru
                theta11 =somru-somru2
                del1i   =somd-somu
                theta11i=somu-somu2
                hpar    =0.
                hpari   =0.
                if(abs(theta11 ).le.tiny(1.)) hpar =del1/theta11
                if(abs(theta11i).le.tiny(1.)) hpari=del1i/theta11i
             else
!             calcul epaisseurs impossible
                del1    =0.
                theta11 =0.
                del1i   =0.
                theta11i=0.
                hpar    =0.
                hpari   =0.
             end if
!
             ii=ideb+idd1*(m1-1)+idd2*(m2-1)
             jj=jdeb+jdd1*(m1-1)+jdd2*(m2-1)
             kk=kdeb+kdd1*(m1-1)+kdd2*(m2-1)
!
!           Mach et Ue parallele a la paroi en delta
             if(iok) then
!             fontiere couche limite
                e=(s(ndel,5)+rpdel*(s(ndel+idm3,5)-s(ndel,5)))/rhoe
                xmue=mu(ndel)+rpdel*(mu(ndel+idm3)-mu(ndel))
!              te=(e-0.5*uemod**2)/cv
                te=temp(ndel)+rpdel*(temp(ndel+idm3)-temp(ndel))
                qq=(usd*usd+utd*utd+und*und)/(rhoe**2)
!              p=gam1*(rhoe*e-.5*rhoe*qq-pinfl)
                p=ps(ndel)+rpdel*(ps(ndel+idm3)-ps(ndel))
                csc=cson(ndel)+rpdel*(cson(ndel+idm3)-cson(ndel))
                sg=sign(1.,qq)
!              rm2=abs(qq*rhoe/(gam*p))
                rm2=abs(qq/csc**2)
                xme=sg*sqrt(rm2)
                uex=usd/rhoe
                reyl=rhoe*sqrt(qq)/xmue
                u=usd
                v=utd
                w=und
             else
!             couche limite non definie. paroi
                n=nc0
                m=n-n0c
                e=s(n,5)/s(n,1)
                qq=(s(n,2)**2+s(n,3)**2+s(n,4)**2)/s(n,1)**2
                p=gam1*(s(n,5)-.5*s(n,1)*qq-pinfl)
                sg=sign(1.,qq*s(n,1)/(gam*p))
                rm2=abs(qq*s(n,1)/(gam*p))
                xme=sg*sqrt(rm2)
                uex=us(m)/s(n,1)
                reyl=s(n,1)*sqrt(qq)/mu(n)
                u=us(m)
                v=ut(m)
                w=un(m)
             endif
!
             write(sor3,'(11(1pe12.4),3i4)') &
                  xcfac,ycfac,zcfac,theta11,hpar,theta11i,hpari, &
                  uex,sqrt(qq),xme,reyl,min(999,ii),min(999,jj),min(999,kk)
!            write(sor3,'(12(1pe12.4),3i4)')
!     &      xcfac,ycfac,zcfac, theta11,hpar, theta11i,hpari,dev,
!     &      uex,sqrt(qq),xme,reyl,min(999,ii),min(999,jj),min(999,kk)
          enddo   !fin de boucle sur les cellules de la bande
       enddo      !fin de boucle sur les bandes
    enddo  !fin de boucle sur les parois
!
    close(sor3)
!
    DEALLOCATE(dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz,vort)

    return
  end subroutine met_intep3
end module mod_met_intep3
