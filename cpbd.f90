      subroutine cpbd( &
                 ncin,nxn,nyn,nzn,ncbd, &
                 sn,vol,v,mut, &
                 bceqt, &
                 rpi,rti,d0x,d0y,d0z,qtx,qty,qtz,x,y,z,omg, &
                 pres,tp,rod,roud,rovd,rowd,roed, &
                 mnr,xnr,ynr,znr,mnc, &
                 tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10,tm11, &
                 tm12,tm13,pression,ztemp,cson)
!
!***********************************************************************
!
!     ACT
!_A    Remplissage des valeurs aux centres d'une rangee de mailles fictives
!_A    tout autour des domaines, par les valeurs aux centres des facettes
!_A    frontieres et les valeurs sur les aretes des maillages. Ces valeurs
!_A    resultent soit d'application des conditions aux limites, soit d'ext-
!_A    trapolations, en particulier en fonction de la cle kexl.
!_A       kexl=1 --->conditions aux limites
!_A       kexl=0 --->extrapolations
!
!     INP
!_I    kexl       : arg int              ; cle traitement frontieres
!_I                                        avant exploitation
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    u          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    rpi        : arg real(ip40      ) ; pres d'arret/pres d'arret etat de
!_I                                        ref utilisateur a imposer
!_I    rti        : arg real(ip40      ) ; temp d'arret/temp d'arret etat de
!_I                                        ref utilisateur a imposer
!_I    d0x        : arg real(ip40      ) ; composante en x d'une
!_I                                        direction de l'ecoulement
!_I    d0y        : arg real(ip40      ) ; composante en y d'une
!_I                                        direction de l'ecoulement
!_I    d0z        : arg real(ip40      ) ; composante en z d'une
!_I                                        direction de l'ecoulement
!_I    qtx        : arg real(ip40      ) ; composante en x d'une
!_I                                        vitesse tangentielle
!_I    qty        : arg real(ip40      ) ; composante en y d'une
!_I                                        vitesse tangentielle
!_I    qtz        : arg real(ip40      ) ; composante en z d'une
!_I                                        vitesse tangentielle
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    omg        : arg real             ; vitesse rotation du repere relatif
!_I    pres       : arg real(ip40      ) ; pression statique
!_I    it         : arg int              ; cycle courant du calcul
!_I    mnr        : arg int (ip44      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule recouvrante
!_I    xnr        : arg real(ip44      ) ; coefficient d'interpolation en x
!_I                                        dans une cellule recouvrante
!_I    ynr        : arg real(ip44      ) ; coefficient d'interpolation en y
!_I                                        dans une cellule recouvrante
!_I    znr        : arg real(ip44      ) ; coefficient d'interpolation en z
!_I                                        dans une cellule recouvrante
!_I    mnc        : arg int (ip43      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule coincidente
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!
!_/    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!_/    cl         : com char(mtb       ) ; type de cond lim a appliquer
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use boundary
      use schemanum
      use sortiefichier
!
!-----------------------------------------------------------------------
!
      character clsave
      real nxn,nyn,nzn
      real mut
!
      dimension tm1(ip40),tm2(ip40),tm3(ip40),tm4(ip40),tm5(ip40), &
                tm6(ip40),tm7(ip40),tm8(ip40),tm9(ip40),tm10(ip40), &
                tm11(ip40),tm12(ip40),tm13(ip40)
      dimension bceqt(ip41,neqt)
      dimension rpi(ip40),rti(ip40),pres(ip40),tp(ip40)
      dimension d0x(ip40),d0y(ip40),d0z(ip40)
      dimension qtx(ip40),qty(ip40),qtz(ip40)
      dimension rod(ip40),roud(ip40),rovd(ip40),rowd(ip40),roed(ip40)
      dimension x(ip21),y(ip21),z(ip21), &
                v(ip11,ip60)
      dimension sn(ip31*ndir)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
      dimension mnr(ip44),xnr(ip44),ynr(ip44),znr(ip44)
      dimension ncin(ip41),mnc(ip43),mut(ip12)
      dimension pression(ip11),ztemp(ip11),cson(ip11),vol(ip11)
      dimension clsave(mtb)
!
      if (kexl.eq.0) then
!
      do no=1,mtbx
       mfb=nba(no)
       clsave(mfb)=cl(mfb)
      enddo
!
      do no=1,mtbx
       mfb=nba(no)
       if (cl(mfb).eq.'idd ') cl(mfb)='rien'
       if (cl(mfb).eq.'idi ') cl(mfb)='rien'
       if (cl(mfb).eq.'iti ') cl(mfb)='rien'
       if (cl(mfb).eq.'gli1') cl(mfb)='rien'
       if (cl(mfb).eq.'glis') cl(mfb)='rien'
       if (cl(mfb).eq.'eikx') cl(mfb)='rien'
       if (cl(mfb).eq.'eikm') cl(mfb)='rien'
       if (cl(mfb).eq.'eikn') cl(mfb)='rien'
       if (cl(mfb).eq.'pari') cl(mfb)='rien'
       if (cl(mfb).eq.'para') cl(mfb)='rien'
       if (cl(mfb).eq.'prec') cl(mfb)='rien'
       if (cl(mfb).eq.'prd ') cl(mfb)='rien'
       if (cl(mfb).eq.'prdp') cl(mfb)='rien'
       if (cl(mfb).eq.'nrd ') cl(mfb)='rien'
       if (cl(mfb).eq.'vrt ') cl(mfb)='rien'
       if (cl(mfb).eq.'axe ') cl(mfb)='rien'
!     condition aux limites lois de paroi
!     lp2 --> lois de paroi standard - parois adiabatiques
!     lp3 --> lois de paroi standard - parois isothermes
!     lp4 --> lois de paroi TBLE 2D  - parois adiabatiques
!     lp5 --> lois de paroi TBLE 3D  - parois adiabatiques
       if (cl(mfb).eq.'lp2 ') cl(mfb)='rien'
       if (cl(mfb).eq.'lp3 ') cl(mfb)='rien'
       if (cl(mfb).eq.'lp4 ') cl(mfb)='rien'
       if (cl(mfb).eq.'lp5 ') cl(mfb)='rien'
       if (cl(mfb).eq.'choc') cl(mfb)='rien'
       if (cl(mfb).eq.'acou') cl(mfb)='rien'
       if (cl(mfb).eq.'debi') cl(mfb)='rien'
      enddo
!
      endif
!
            call rbc( &
                 ncin,nxn,nyn,nzn,ncbd, &
                 sn,vol,v,mut, &
                 bceqt, &
                 rpi,rti,d0x,d0y,d0z,qtx,qty,qtz,x,y,z,omg, &
                 pres,tp,rod,roud,rovd,rowd,roed, &
                 numt, &
                 mnr,xnr,ynr,znr,mnc, &
                 tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10,tm11, &
                 tm12,tm13,pression,ztemp,cson)
!
      do no=1,mtbx
       mfb=nba(no)
       cl(mfb)=clsave(mfb)
      enddo
!
      return
      end
