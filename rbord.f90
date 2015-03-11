module mod_rbord
implicit none
contains
      subroutine rbord( &
                 kpst,img, &
                 mnc,ncin,mnr,xnr,ynr,znr,nxn,nyn,nzn,ncbd, &
                 sn,vol,v,dt, &
                 bceqt, &
                 rpi,rti,d0x,d0y,d0z,x,y,z, &
                 pres,tp,rod,roud,rovd,rowd,roed, &
                 icyc, &
                 tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10,tm11, &
                 tm12,tm13, &
                 pression,ztemp,cson)
!
!***********************************************************************
!
!     ACT
!_A    Application des conditions aux limites sur chacune des frontieres
!_A    limites. Loi d'etat des gaz raides. le fluide est de l'eau liquide.
!
!     INP
!_I    kpst       : arg int              ; cle traitement de la condition limite
!_I                                        'rien' en cas de post-traitement
!_I    mnc        : arg int (ip43      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule coincidente
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    mnr        : arg int (ip44      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule recouvrante
!_I    xnr        : arg real(ip44      ) ; coefficient d'interpolation en x
!_I                                        dans une cellule recouvrante
!_I    ynr        : arg real(ip44      ) ; coefficient d'interpolation en y
!_I                                        dans une cellule recouvrante
!_I    znr        : arg real(ip44      ) ; coefficient d'interpolation en z
!_I                                        dans une cellule recouvrante
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
!_I    pres       : arg real(ip40      ) ; pression statique
!_I    tp         : arg real(ip40      ) ; temperature paroi a imposer
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    it         : arg int              ; cycle courant du calcul
!_I    cl         : com char(mtb       ) ; type de cond lim a appliquer
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    omg        : com real             ; vitesse rotation du repere relatif
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    nnn        : com int (lt        ) ; nombre de noeuds du dom (dont fic.)
!_I    npfb       : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes facettes
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use sortiefichier
      use maillage
      use boundary
      use chainecarac
      use schemanum
      use definition
      use modeleturb
use mod_utprd
use mod_utpari
use mod_clgli2
use mod_clprd
use mod_clpara
use mod_clidd2_prcd
use mod_clpari
use mod_clnrd_prcd
use mod_clidd
use mod_clchoc
use mod_lecture_acou
use mod_utvrt
use mod_clextr
use mod_rbvr
use mod_clglis
use mod_clsym
use mod_clnrd
use mod_rbvc
use mod_utnrd
use mod_cldebit
use mod_rbve
use mod_utdebit
use mod_utidd
use mod_clidi0
use mod_utidi
use mod_clvrt
use mod_clprd_prcd
use mod_cldebit_prcd
implicit none
integer :: kpst
integer :: img
integer :: mnc
integer :: ncin
integer :: mnr
double precision :: xnr
double precision :: ynr
double precision :: znr
integer :: ncbd
double precision :: sn
double precision :: vol
double precision :: v
double precision :: dt
double precision :: bceqt
double precision :: rpi
double precision :: rti
double precision :: d0x
double precision :: d0y
double precision :: d0z
double precision :: x
double precision :: y
double precision :: z
double precision :: pres
double precision :: tp
double precision :: rod
double precision :: roud
double precision :: rovd
double precision :: rowd
double precision :: roed
integer :: icyc
double precision :: tm1
double precision :: tm2
double precision :: tm3
double precision :: tm4
double precision :: tm5
double precision :: tm6
double precision :: tm7
double precision :: tm8
double precision :: tm9
double precision :: tm10
double precision :: tm11
double precision :: tm12
double precision :: tm13
double precision :: pression
double precision :: ztemp
double precision :: cson
integer :: l
integer :: lgsnlt
integer :: lm
integer :: mfb
integer :: mfbm
integer :: no
integer :: npsn
!
!-----------------------------------------------------------------------
!
      real nxn,nyn,nzn
      dimension tm1(ip40),tm2(ip40),tm3(ip40),tm4(ip40),tm5(ip40), &
                tm6(ip40),tm7(ip40),tm8(ip40),tm9(ip40),tm10(ip40), &
                tm11(ip40),tm12(ip40),tm13(ip40)
      dimension bceqt(ip41,neqt)
      dimension rpi(ip40),rti(ip40),pres(ip40),tp(ip40)
      dimension d0x(ip40),d0y(ip40),d0z(ip40)
      dimension rod(ip40),roud(ip40),rovd(ip40),rowd(ip40),roed(ip40)
      dimension ncin(ip41),mnc(ip43)
      dimension x(ip21),y(ip21),z(ip21)
      dimension v(ip11,ip60)
      dimension sn(ip31*ndir)
      dimension vol(ip11),pression(ip11),ztemp(ip11),cson(ip11),dt(ip11)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
      dimension xnr(ip44),ynr(ip44),znr(ip44),mnr(ip44)
!
      nbd=1
!
      do no=1,mtbx
!
      mfb   =nba(no)
      mfbm  =mfb+(img-1)*mtb
      lbd(1)=mfbm
      l     =ndlb(mfb)
      lm    =l+(img-1)*lz
!
!-----condition d' injection , direction de la vitesse imposee
      if (cl(mfb)(1:3).eq.'idd') then
            call utidd ( &
                 bceqt, &
                 mfbm,rpi,rti,d0x,d0y,d0z,ncbd, &
                 mmb,mpb)
            call rbve(v,ncin,ncbd)
            if(kprec.eq.0) then
             call clidd( &
                 mfbm,lm,rpi,rti,d0x,d0y,d0z, &
                 nxn,nyn,nzn,ncbd,v, &
                 y,z, &
                 mmb,mpb,mpn, &
                 tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9, &
                 tm10,tm11,pression,ztemp,cson)
            elseif(kprec.ge.1) then
             call clidd2_prcd( &
                 mfbm,lm,rpi,rti,d0x,d0y,d0z, &
                 nxn,nyn,nzn,ncbd,v, &
                 y,z, &
                 mmb,mpb,mpn, &
                 tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9, &
                 tm10,tm11,tm12,pression,ztemp,cson)
            endif
!
!-----condition d' injection , (direction de la vitesse)/dt=0
      else if (cl(mfb)(1:3).eq.'idi') then
            call utidi ( &
                 bceqt, &
                 mfbm,rpi,rti,ncbd, &
                 mmb,mpb)
            call rbve(v,ncin,ncbd)
!
            call clidi0( &
                 mfbm,lm,rpi,rti, &
                 nxn,nyn,nzn,ncbd,v, &
                 y,z, &
                 mmb,mpb,mpn, &
                 tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10,tm11, &
                 tm12,tm13,pression,ztemp,cson)
!
!-----condition de glissement
      else if (cl(mfb).eq.'gli1') then
             call clgli2( &
                 ncbd,ncin, &
                 mfbm,mmb,mpb,mpn, &
                 nxn,nyn,nzn, &
                 v,pression,ztemp,cson)
!
!-----condition de glissement complique!
      else if (cl(mfb).eq.'glis') then
            call rbve(v,ncin,ncbd)
!
      npsn  =ndir*npfb(lm)+1
      lgsnlt=nnn(lm)
!
            call clglis( &
                 mfbm,lm,indfl(mfb), &
                 ncin,ncbd, &
                 x,y,z, &
                 sn(npsn),lgsnlt,vol, &
                 v,pression,ztemp,cson)
!
!-----condition de paroi isotherme
      else if (cl(mfb).eq.'pari') then
            call utpari( &
                 bceqt, &
                 mfbm,tp, &
                 mmb,mpb)
            call clpari( &
                 mfbm,tp, &
                 ncbd,v, &
                 mmb,mpb,ncin, &
                 pression,ztemp,cson)
!
!-----condition de paroi adiabatique
      else if (cl(mfb).eq.'para') then
            call clpara( &
                 mfbm, &
                 ncbd,v, &
                 mmb,mpb,ncin, &
                 pression,ztemp,cson)
!
!-----condition de pression , pression imposee(relations de compatibilite)
      else if (cl(mfb)(1:3).eq.'prd') then
            call utprd( &
                 bceqt, &
                 mfbm,pres,icyc, &
                 mmb,mpb)
            call rbve(v,ncin,ncbd)
!
            if(kprec.eq.0) then
               call clprd( &
                 mfbm,pres, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,lm, &
                 pression,ztemp,cson)
!
            elseif(kprec.ge.1) then
               call clprd_prcd( &
                 mfbm,pres, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,lm, &
                 pression,ztemp,cson)
            endif
!
!-----condition de debit, debit impose en sortie
      else if (cl(mfb)(1:4).eq.'debi') then
            call rbve(v,ncin,ncbd)
!
            npsn  =ndir*npfb(lm)+1
            lgsnlt=nnn(lm)
!
            call utdebit( &
                 mfbm,bceqt,pres, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,lm, &
                 sn(npsn),lgsnlt)
!
            if(kprec.eq.0) then
             call cldebit( &
                 mfbm,pres, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,lm, &
                 pression,ztemp,cson)
            elseif(kprec.ge.1) then
             call cldebit_prcd( &
                 mfbm,pres, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,lm, &
                 pression,ztemp,cson)
            endif
!
!-----condition d'extrapolation
      else if (cl(mfb).eq.'extr') then
            call clextr( &
                 mfbm, &
                 ncbd,v, &
                 mmb,mpb,ncin, &
                 pression,ztemp,cson)
!
!-----generation d'un choc oblique
      else if (cl(mfb).eq.'choc') then
            call clchoc( &
                 mfbm,mmb,mpb, &
                 ncbd,v)
!
!-----condition de non reflexion
      else if (cl(mfb)(1:3).eq.'nrd') then
            call utnrd( &
                 bceqt, &
                 mfbm,lm,rod,roud,rovd,rowd,roed, &
                 ncbd, &
                 y,z)
            call rbve(v,ncin,ncbd)
!
            if(kprec.eq.0) then
             call clnrd( &
                 mfbm,rod,roud,rovd,rowd,roed, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,lm, &
                 pression,ztemp,cson)
!
            elseif(kprec.ge.1) then
             call clnrd_prcd( &
                 mfbm,rod,roud,rovd,rowd,roed, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,lm, &
                 pression,ztemp,cson)
            endif
!
!-----condition d'atmosphere infinie pour maillage "10 cordes" en transsonique
      else if (cl(mfb)(1:3).eq.'vrt') then
            call utvrt( &
                 bceqt, &
                 mfbm,lm,rod,roud,rovd,rowd,roed, &
                 ncbd, &
                 y,z, &
                 mmb,mpb)
            call rbve(v,ncin,ncbd)
            call clvrt( &
                 mfbm,rod,roud,rovd,rowd,roed, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn, &
                 l,x,y,z, &
                 pression,ztemp,cson)
!
!-----condition de symetrie par rapport aux facettes frontieres
      else if (cl(mfb)(1:3).eq.'sym') then
            call clsym( &
                 mfbm, &
                 nxn,nyn,nzn,ncin,ncbd,v, &
                 mmb,mpb,mpn)
!
!-----raccord recouvrant
      else if (cl(mfb)(1:3).eq.'rec') then
            call rbvr( &
                 ncbd,mnr,xnr,ynr,znr, &
                 v, &
                 ncin)
!
!-----raccord non coincident
      else if (cl(mfb)(1:3).eq.'rnc') then
            call rbvr( &
                 ncbd,mnr,xnr,ynr,znr, &
                 v, &
                 ncin)
!
!-----raccord coincident
      else if (cl(mfb)(1:2).eq.'rc') then
            call rbvc(v,ncbd,ncin,mnc)
!
!-----affectation de valeurs extrapolees en vue de post-traitement
      else if (((cl(mfb).eq.'rien').or.(cl(mfb).eq.'axe ')).and. &
               (kpst.eq.1)) then
            call rbve(v,ncin,ncbd)
!
!-----aucun traitement n'est effectue
      else if (((cl(mfb).eq.'rien').or.(cl(mfb).eq.'axe ')).and. &
               (kpst.ne.1)) then
!
!-----lois de paroi
      else if(cl(mfb)(1:2).eq.'lp') then
!       constantes lois de paroi
        vkar=0.41
        cllog=5.25
        yp0=11.13
        if(cl(mfb)(1:3).eq.'lp2') then
!         cle lois de paroi
          lparoi=1
!         parois adiabatiques
          call clpara( &
                 mfbm, &
                 ncbd,v, &
                 mmb,mpb,ncin, &
                 pression,ztemp,cson)
        elseif(cl(mfb)(1:3).eq.'lp3') then
!         cle lois de paroi
          lparoi=1
!         parois isothermes
          call utpari( &
                 bceqt, &
                 mfbm,tp, &
                 mmb,mpb)
            call clpari( &
                 mfbm,tp, &
                 ncbd,v, &
                 mmb,mpb,ncin, &
                 pression,ztemp,cson)
        elseif((cl(mfb)(1:3).eq.'lp4').or.(cl(mfb)(1:3).eq.'lp5')) then
!         cle lois de paroi
          lparoi=2
!         parois adiabatiques
          call clpara( &
                 mfbm, &
                 ncbd,v, &
                 mmb,mpb,ncin, &
                 pression,ztemp,cson)
        endif
!
        elseif (cl(mfb)(1:4).eq.'acou') then
           if(icyc.eq.1) then
            call lecture_acou ( &
                 1,v)
           endif
      else
      if (kimp.ge.1) then
      write(imp,'(a,a)') 'condition limite inconnue : ',cl(mfb)
      endif
      stop
      end if
!
      enddo
!
      return
      end subroutine
end module
