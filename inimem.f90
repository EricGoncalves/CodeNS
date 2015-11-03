module mod_inimem
  implicit none
contains
  subroutine inimem( &
       dt,v,mu,mut, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       sn, &
       vol, &
       ptdual,vdual,vdual1,vdual2, &
       cvi,cvj,cvk, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
       pression,ztemp,cson, &
       tnte1,tnte3,tnte4, &
       ncyc, &
       tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
       comment)
!
!***********************************************************************
!
!     ACT
!_A    Mise a vide de la variable de commentaire.
!_A    Mise a zero de l'encombrement effectif des tableaux.
!_A    Initialisation des valeurs par defaut.
!_A    Mise a zero du numero de cycle courant.
!_A    Initialisation des tableaux pour vectorisation.
!
!_O    u          : arg real(ip11,ip60 ) ; variables a l'instant n
!_O    dt         : arg real(ip11      ) ; pas de temps
!_O    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!_O    mut        : arg real(ip12      ) ; viscosite turbulente
!_O    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_O                                        contraintes visqueuses
!_O    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_O                                        contraintes visqueuses
!_O    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_O                                        contraintes visqueuses
!_O    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_O                                        contraintes visqueuses
!_O    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_O                                        contraintes visqueuses
!_O    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_O                                        contraintes visqueuses
!_O    qcx        : arg real(ip12      ) ; composante en x du flux de chaleur
!_O    qcy        : arg real(ip12      ) ; composante en y du flux de chaleur
!_O    qcz        : arg real(ip12      ) ; composante en z du flux de chaleur
!_O    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_O                                        norme egale a la surface de celle-ci
!_O    vol        : arg real(ip11      ) ; volume d'une cellule
!_O    ncyc       : arg int              ; cycle courant de l'execution
!_O    equat      : com char             ; type d'equations modelisant l'ecoule-
!_O                                        ment
!_O    lzx        : com int              ; nbr total de domaines
!_O    mtbx       : com int              ; nbr total de frontieres
!_O    mtnx       : com int              ; nbr total de frontieres a normales stockes
!_O    mtcx       : com int              ; nbr total de frontieres coincidentes
!_O    mtrx       : com int              ; nbr total de frontieres recouvertes
!_O    mtax       : com int              ; nbr total de frontieres autres
!_O    ndimubx    : com int              ; nbr de cellules du plus grd domaine
!_O                                        (pts fictifs inclus)
!_O    ndimctbx   : com int              ; nbr de cellules de tts les domaines
!_O                                        (pts fictifs inclus)
!_O    ndimntbx   : com int              ; nbr de noeuds de tts les domaines
!_O                                        (pts fictifs inclus)
!_O    mdimubx    : com int              ; nbr de pts de la plus grde front
!_O    mdimtbx    : com int              ; nbr de pts de ttes les front
!_O    mdimtnx    : com int              ; nbr de pts de ttes les front a normales stockees
!_O    mdimtcx    : com int              ; nbr de pts de ttes les front coincidentes
!_O    mdimtrx    : com int              ; nbr de pts de ttes les front recouvertes
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use kcle
    use schemanum
    use proprieteflu
    use mod_def
    use constantes 
    use mod_inivec
    implicit none
    integer          :: icmt,ncyc
    double precision ::    cmui1,   cmui2,   cmuj1,   cmuj2,   cmuk1
    double precision ::    cmuk2,    cson,     cvi,     cvj,     cvk
    double precision ::       dt,      mu,     mut,pression,  ptdual
    double precision ::      qcx,     qcy,     qcz,      sn,     tn1
    double precision ::     tn10,     tn2,     tn3,     tn4,     tn5
    double precision ::      tn6,     tn7,     tn8,     tn9,   tnte1
    double precision ::    tnte3,   tnte4,    toxx,    toxy,    toxz
    double precision ::     toyy,    toyz,    tozz,       v,   vdual
    double precision ::   vdual1,  vdual2,     vol,   ztemp
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
!
    dimension v(ip11,ip60)
    dimension dt(ip11),vol(ip11),pression(ip11),ztemp(ip11),cson(ip11)
    dimension mu(ip12),mut(ip12)
    dimension toxx(ip12),toxy(ip12),toxz(ip12),toyy(ip12),toyz(ip12), &
         tozz(ip12),qcx(ip12),qcy(ip12),qcz(ip12)
    dimension sn(ip31*ndir)
    dimension tn1(ip00),tn2(ip00),tn3(ip00),tn4(ip00),tn5(ip00),tn6(ip00), &
         tn7(ip00),tn8 (ip00),tn9(ip00),tn10(ip00)
    dimension tnte1(ip11,ip60),tnte3(ip11,ip60),tnte4(ip11,ip60)
    dimension vdual(ip11,ip60),vdual1(ip11,ip60),vdual2(ip11,ip60),ptdual(ip11,ip60)
    dimension cvi(ip21),cvj(ip21),cvk(ip21), &
         cmui1(ip21),cmui2(ip21),cmuj1(ip21),cmuj2(ip21), &
         cmuk1(ip21),cmuk2(ip21)
!
!     mise a vide de la variable de commentaire
!
    do icmt=1,32
       comment(icmt:icmt)=' '
    enddo
!
!     dimensionnement du code
!
    lzx      =0
    lgx      =0
    mtbx     =0
    mtrx     =0
    mtcx     =0
    mtnx     =0
    mtax     =0
    ndimubx  =0
    ndimctbx=0
    ndimntbx =0
    mdimubx  =0
    mdimtbx  =0
    mdimtnx  =0
    mdimtcx  =0
    mdimtrx  =0
!
    klzx      =0
    kmtbx     =0
    kmtrx     =0
    kmtcx     =0
    kmtnx     =0
    kmtax     =0
    kndimubx  =0
    kndimctbx =0
    kndimntbx =0
    kmdimubx  =0
    kmdimtbx  =0
    kmdimtnx  =0
    kmdimtcx  =0
    kmdimtrx  =0

!     initialisation des constantes

    intmx=huge(1)
    reelmx=huge(1.)
    reelmn=tiny(1.)
    pis2=atan2(1.D0,0.D0)
    raddeg=90.D0/pis2
    degrad=pis2/90.D0
!
!     valeurs par defaut
!
    call def
!
!     nombre de cycles precedents en debut d'execution
    ncyc=0
    tol=1.e-08
    tolke=0.
    niter=1
    nitur=1
    kdualns=0
    kdualto=0
    resno1=1.
    resite=1.
    ischema=1
    muscl=0
    ilim=0
    xk=1.
    kprec=0
    cte=0.5
    kvisq=0
    ql=0.
    pinfl=1.
!
!     initialisation des tableaux
    call inivec( &
         dt,v,mu,mut, &
         toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
         sn, &
         vol, &
         ptdual,vdual,vdual1,vdual2, &
         cvi,cvj,cvk, &
         cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
         pression,ztemp,cson, &
         tnte1,tnte3,tnte4, &
         tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10)
!
    return
  end subroutine inimem
end module mod_inimem
