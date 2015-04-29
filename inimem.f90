module mod_inimem
  implicit none
contains
  subroutine inimem(ncyc)
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
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  comment
!
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
    pis2=acos(-1.)/2.
    raddeg=90./pis2
    degrad=pis2/90.

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
    return
  end subroutine inimem
end module mod_inimem
