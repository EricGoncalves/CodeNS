module mod_b2_dpdim
  implicit none
contains
  subroutine b2_dpdim
!
!***********************************************************************
!
!     ACT
!_A    Sorties dans le fichier fimp pour la commande 'dpdim'.
!
!_I    kdimg      : arg int              ; dim, cle tab multigrille
!_I    kdimv      : arg int              ; dim, cle tab visqueux
!_I    kdimk      : arg int              ; dim, cle tab k-eps
!_I    equat      : com char             ; type d'equations modelisant l'ecoulement
!_I    lzx        : com int              ; nbr total de domaines
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    ndimubx    : com int              ; nbr de cellules du plus grd domaine
!_I                                        (pts fictifs inclus)
!_I    ndimctbx   : com int              ; nbr de cellules de tts les domaines
!_I                                        (pts fictifs inclus)
!_I    ndimntbx   : com int              ; nbr de noeuds de tts les domaines
!_I                                        (pts fictifs inclus)
!_I    mdimubx    : com int              ; nbr de pts de la plus grde front
!_I    mdimtbx    : com int              ; nbr de pts de ttes les front
!_I    mdimtnx    : com int              ; nbr de pts de ttes les front
!_I                                        a normales stockees
!_I    mdimtcx    : com int              ; nbr de pts de ttes les front coincidentes
!_I    mdimtrx    : com int              ; nbr de pts de ttes les front recouvertes
!_I    kmf        : com int (lt        ) ; cle phase implicite
!_I    imp        : com int              ; unite logiq, sorties de controle
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use schemanum
    use maillage
    use chainecarac
    implicit none
    integer          ::      kl,      l, mdimtb, mdimtc, mdimtn
    integer          ::  mdimtr,ndimctb,ndimctc,ndimctk,ndimctv
    integer          :: ndimntb
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
!
    ndimctb=nint((1.+kdimg*ccg)*ndimctf)
    ndimctv=kdimv*(ndimctb-1)+1
    ndimctk=kdimk*(ndimctb-1)+1
    ndimctc=nint(kdimg*ccg*ndimctf)
    ndimntb=nint((1.+kdimg*cng)*ndimnts+ndimntu)
    mdimtb =nint((1.+kdimg*cfg)*mdimtbf)
    mdimtn =nint((1.+kdimg*cfg)*mdimtnf)
    mdimtc =nint((1.+kdimg*cfg)*mdimtcf)
    mdimtr =nint((1.+kdimg*cfg)*mdimtrf)
!
    kl=0
    do l=1,lzx
       kl=max(kl,kmf(l))
    enddo
!
    form='(''Verification du dimensionnement du code:'',/,' &
         //'1x,39(1h-),                             //,' &
         //'1x,''code     '',  11x,''calcul    '',7x/,' &
         //'1x,''----     '',  11x,''------    '',7x//,' &
         //'1x,''lz     = '',i7,4x,''lzx     = '',i7//,' &
         //'1x,''ndimub = '',i7,4x,''ndimubx = '',i7//,' &
         //'1x,''ndimctb= '',i7,4x,''ndimctbx= '',i7//,' &
         //'1x,''ndimntb= '',i7,4x,''ndimntbx= '',i7//,' &
         //'1x,''nvar   = '',i7,4x,''equat   = '',2x,a//,' &
         //'1x,''kdimv  = '',i7,//,' &
         //'1x,''kdimk  = '',i7,//,' &
         //'1x,''kdimg  = '',i7,4x,''lgx     = '',i7/)'
    write(imp,form)lz,lzx,ndimub,ndimubx,ndimctb,ndimctbx, &
         ndimntb,ndimntbx, &
         nvar,equat,kdimv,kdimk,kdimg,lgx
    form='(' &
         //'1x,''mtb     = '',i7,4x,''mtbx    = '',i7//,' &
         //'1x,''mdimub  = '',i7,4x,''mdimubx = '',i7//,' &
         //'1x,''mdimtb  = '',i7,4x,''mdimtbx = '',i7//,' &
         //'1x,''mdimtn  = '',i7,4x,''mdimtnx = '',i7//,' &
         //'1x,''mdimtc  = '',i7,4x,''mdimtcx = '',i7//,' &
         //'1x,''mdimtr  = '',i7,4x,''mdimtrx = '',i7//)'
    write(imp,form)mtb,mtbx, &
         mdimub,mdimubx,mdimtb,mdimtbx,mdimtn,mdimtnx, &
         mdimtc,mdimtcx,mdimtr,mdimtrx
!
    return
  end subroutine b2_dpdim
end module mod_b2_dpdim
