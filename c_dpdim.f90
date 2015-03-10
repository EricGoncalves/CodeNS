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

implicit none
integer :: imot
integer :: nmot
integer :: kl
integer :: l
integer :: mdimtb
integer :: mdimtc
integer :: mdimtn
integer :: mdimtr
integer :: ndimctb
integer :: ndimctc
integer :: ndimctk
integer :: ndimctv
integer :: ndimntb
!
!-----------------------------------------------------------------------
!
      character(len=32) ::  mot(nmx)
      dimension imot(nmx)
!
      ndimctb=int((1.+kdimg*ccg)*ndimctf)
      ndimctv=kdimv*(ndimctb-1)+1
      ndimctk=kdimk*(ndimctb-1)+1
      ndimctc=int(kdimg*ccg*ndimctf)
      ndimntb=int((1.+kdimg*cng)*ndimnts+ndimntu)
      mdimtb =int((1.+kdimg*cfg)*mdimtbf)
      mdimtn =int((1.+kdimg*cfg)*mdimtnf)
      mdimtc =int((1.+kdimg*cfg)*mdimtcf)
      mdimtr =int((1.+kdimg*cfg)*mdimtrf)
!
      kl=0
      do l=1,lzx
       kl=max(kl,kmf(l))
      enddo
!
      if (kimp.ge.2) then
            call b2_dpdim
      endif
!
      if (lt     .lt.lzx     ) stop 'dimensionnement incorrecte'
      if (ndimub .lt.ndimubx ) stop 'dimensionnement incorrecte'
      if (ndimctb.lt.ndimctbx) stop 'dimensionnement incorrecte'
      if (ndimntb.lt.ndimntbx) stop 'dimensionnement incorrecte'
      if((kdimg.eq.0).and.(lgx.gt.1)) &
                               stop 'dimensionnement incorrecte'
      if((kdimv.eq.0).and.(equat(1:2).eq.'ns')) &
                               stop 'dimensionnement incorrecte'
      if((kdimk.eq.0).and.(equat(1:2).eq.'ke')) &
                               stop 'dimensionnement incorrecte'
      if (mtb    .lt.mtbx    ) stop 'dimensionnement incorrecte'
      if (mdimub .lt.mdimubx ) stop 'dimensionnement incorrecte'
      if (mdimtb .lt.mdimtbx ) stop 'dimensionnement incorrecte'
      if (mdimtn .lt.mdimtnx ) stop 'dimensionnement incorrecte'
      if (mdimtc .lt.mdimtcx ) stop 'dimensionnement incorrecte'
      if (mdimtr .lt.mdimtrx ) stop 'dimensionnement incorrecte'
!
      return
      end
end module
