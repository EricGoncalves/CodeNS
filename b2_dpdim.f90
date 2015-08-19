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
    use mod_mpi
    implicit none
    integer          ::      kl,      l, mdimtb, mdimtc, mdimtn
    integer          ::  mdimtr,ndimctb,ndimctc,ndimctk,ndimctv
    integer          :: ndimntb,buff(6)
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
    call SUM_MPI(lzx,BUFF(1))
    call MAX_MPI(ndimubx,BUFF(2))
    call MAX_MPI(ndimctbx,BUFF(3))
    call MAX_MPI(ndimntbx,BUFF(4))

    form='(''Verification du dimensionnement du code:'',/,' &
         //'1x,39(1h-),                             //,' &
         //'1x,''code     '',  11x,''calcul           '',7x/,' &
         //'1x,''----     '',  11x,''------           '',7x//,' &
         //'1x,''lz     = '',i7,4x,''sum of lzx     = '',i7//,' &
         //'1x,''ndimub = '',i7,4x,''max of ndimubx = '',i7//,' &
         //'1x,''ndimctb= '',i7,4x,''max of ndimctbx= '',i7//,' &
         //'1x,''ndimntb= '',i7,4x,''max of ndimntbx= '',i7//,' &
         //'1x,''nvar   = '',i7,4x,''equat          = '',2x,a//,' &
         //'1x,''kdimv  = '',i7,//,' &
         //'1x,''kdimk  = '',i7,//,' &
         //'1x,''kdimg  = '',i7,4x,''lgx            = '',i7/)'
    if (rank==0)write(imp,form)lz,BUFF(1),ndimub,BUFF(2),ndimctb,BUFF(3), &
         ndimntb,BUFF(4), &
         nvar,equat,kdimv,kdimk,kdimg,lgx

    call MAX_MPI(mtbx,BUFF(1))
    call MAX_MPI(mdimubx,BUFF(2))
    call MAX_MPI(mdimtbx,BUFF(3))
    call MAX_MPI(mdimtnx,BUFF(4))
    call MAX_MPI(mdimtcx,BUFF(5))
    call MAX_MPI(mdimtrx,BUFF(6))

    form='(' &
         //'1x,''mtb     = '',i7,4x,''max of mtbx    = '',i7//,' &
         //'1x,''mdimub  = '',i7,4x,''max of mdimubx = '',i7//,' &
         //'1x,''mdimtb  = '',i7,4x,''max of mdimtbx = '',i7//,' &
         //'1x,''mdimtn  = '',i7,4x,''max of mdimtnx = '',i7//,' &
         //'1x,''mdimtc  = '',i7,4x,''max of mdimtcx = '',i7//,' &
         //'1x,''mdimtr  = '',i7,4x,''max of mdimtrx = '',i7//)'
    if (rank==0)write(imp,form)mtb,BUFF(1), &
         mdimub,BUFF(2),mdimtb,BUFF(3),mdimtn,BUFF(4), &
         mdimtc,BUFF(5),mdimtr,BUFF(6)
!
    return
  end subroutine b2_dpdim
end module mod_b2_dpdim
