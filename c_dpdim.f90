module mod_c_dpdim
  implicit none
contains
  subroutine c_dpdim(mot,imot,nmot)
!
!***********************************************************************
!
!     ACT
!_A    Realisation de l'action dpdim.
!
!_I    equat      : com char             ; type d'equations modelisant l'ecoule-
!_I                                        ment
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
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
!_I    mdimtcx    : com int              ; nbr de pts de ttes les front
!_I                                        coincidentes
!_I    mdimtrx    : com int              ; nbr de pts de ttes les front
!_I                                        recouvertes
!_I    kmf       : com int (lt        ) ; cle phase implicite
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use chainecarac
    use schemanum
    use boundary
    use sortiefichier
    use mod_b2_dpdim
    use mod_mpi
    implicit none
    integer          :: imot(nmx),       kl,        l,   mdimtb,   mdimtc
    integer          ::    mdimtn,   mdimtr,  ndimctb,  ndimctc,  ndimctk
    integer          ::   ndimctv,  ndimntb,     nmot
!
!-----------------------------------------------------------------------
!
    character(len=32) ::  mot(nmx)
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
    if (kimp.ge.2) then
       call start_keep_order
       write(imp,*) "rank :",rank
       call b2_dpdim
       call end_keep_order
    endif
!
    if (lt     .lt.lzx     ) stop 'dimensionnement incorrect lt'
    if (ndimub .lt.ndimubx ) stop 'dimensionnement incorrect ndimub'
    if (ndimctb.lt.ndimctbx) stop 'dimensionnement incorrect ndimctb'
    if (ndimntb.lt.ndimntbx) stop 'dimensionnement incorrect ndimntb'
    if((kdimg.eq.0).and.(lgx.gt.1)) &
         stop 'dimensionnement incorrect kdimg'
    if((kdimv.eq.0).and.(equat(1:2).eq.'ns')) &
         stop 'dimensionnement incorrect kdimv'
    if((kdimk.eq.0).and.(equat(1:2).eq.'ke')) &
         stop 'dimensionnement incorrect kdimk'
    if (mtb    .lt.mtbx    ) stop 'dimensionnement incorrect mtb'
    if (mdimub .lt.mdimubx ) stop 'dimensionnement incorrect mdimub'
    if (mdimtb .lt.mdimtbx ) stop 'dimensionnement incorrect mdimtb'
    if (mdimtn .lt.mdimtnx ) stop 'dimensionnement incorrect mdimtn'
    if (mdimtc .lt.mdimtcx ) stop 'dimensionnement incorrect mdimtc'
    if (mdimtr .lt.mdimtrx ) stop 'dimensionnement incorrect mdimtr'
!
    return
  end subroutine c_dpdim
end module mod_c_dpdim
