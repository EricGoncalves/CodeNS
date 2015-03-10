module mod_crbds
implicit none
contains
      subroutine crbds( &
                 mfbe,kini,l, &
                 imin,imax,jmin,jmax,kmin,kmax, &
                 indmf, &
                 ncbd)
!
!***********************************************************************
!
!     ACT
!_A    Creation d'une frontiere essentiellement par le reperage des
!_A    points de la frontiere dans les tableaux domaines.
!
!     INP
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    l          : arg int              ; numero de domaine
!_I    kini       : arg int              ; cle creation de frontiere
!_I    imin       : arg int              ; indice min en i
!_I    imax       : arg int              ; indice max en i
!_I    jmin       : arg int              ; indice min en j
!_I    jmax       : arg int              ; indice max en j
!_I    kmin       : arg int              ; indice min en k
!_I    kmax       : arg int              ; indice max en k
!_I    indmf      : arg char             ; type de plan de la frontiere
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    kfi        : com int              ; unite logiq, indices pts frontieres
!
!     OUT
!_O    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_O                                        cellule frontiere fictive
!_O    indfl      : com char(mtb       ) ; type de plan de la frontiere
!_O    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_O    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_O                                        en fct du numero externe
!_O    iminb      : com int (mtt       ) ; indice min en i d'une frontiere
!_O    imaxb      : com int (mtt       ) ; indice max en i d'une frontiere
!_O    jminb      : com int (mtt       ) ; indice min en j d'une frontiere
!_O    jmaxb      : com int (mtt       ) ; indice max en j d'une frontiere
!_O    kminb      : com int (mtt       ) ; indice min en k d'une frontiere
!_O    kmaxb      : com int (mtt       ) ; indice max en k d'une frontiere
!
!     I/O
!_/    mtbx       : com int              ; nbr total de frontieres
!_/    mdimubx    : com int              ; nbr de pts de la plus grde front
!_/    mdimtbx    : com int              ; nbr de pts de ttes les front
!_/    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_/    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_/                                        dans tableaux de base des front.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use kcle
      use chainecarac
      use sortiefichier
      use boundary
use mod_initis
implicit none
integer :: mfbe
integer :: kini
integer :: l
integer :: imin
integer :: imax
integer :: jmin
integer :: jmax
integer :: kmin
integer :: kmax
integer :: ncbd
integer :: img
integer :: imgi
integer :: imgj
integer :: imgk
integer :: lm
integer :: m0
integer :: mfbi
integer :: mfbim
integer :: mt
!
!-----------------------------------------------------------------------
!
      character(len=2 ) :: indmf
      dimension ncbd(ip41)
!
      mtbx=mtbx+1
      kmtbx=2
!
      mfbi=mtbx
      nfei(mfbe)=mfbi
      ndlb(mfbi)=l
      indfl(mfbi)=indmf
!
      do img=1,lgx
!
      lm=l+(img-1)*lz
      mfbim=mfbi+(img-1)*mtb
!
      imgi=img
      imgj=img
      imgk=img
      if (equat(3:5).eq.'2di') imgi = 1
      if (equat(3:5).eq.'2dj') imgj = 1
      if (equat(3:5).eq.'2dk') imgk = 1
      if (equat(3:5).eq.'2xk') imgk = 1
!
      iminb(mfbim)=(imin-ii1(lm))/2**(imgi-1)+ii1(lm)
      imaxb(mfbim)=(imax-ii1(lm))/2**(imgi-1)+ii1(lm)
      jminb(mfbim)=(jmin-jj1(lm))/2**(imgj-1)+jj1(lm)
      jmaxb(mfbim)=(jmax-jj1(lm))/2**(imgj-1)+jj1(lm)
      kminb(mfbim)=(kmin-kk1(lm))/2**(imgk-1)+kk1(lm)
      kmaxb(mfbim)=(kmax-kk1(lm))/2**(imgk-1)+kk1(lm)
!
      mpb(mfbim)=mdimtbx
      m0=mpb(mfbim)
!
!     remplissage des tableaux  ncbd, mmb
!
      if((kini.eq.1).or.(img.gt.1)) then
            call initis( &
                 lm, &
                 iminb(mfbim),imaxb(mfbim), &
                 jminb(mfbim),jmaxb(mfbim), &
                 kminb(mfbim),kmaxb(mfbim), &
                 indmf,ncbd, &
                 mt,m0)
      elseif(kini.eq.0) then
!            call readfi( &
!                 kfi,ncbd, &
!                 mt,m0)
      endif
!
      mmb(mfbim)=mt
      mdimubx=max(mdimubx,mmb(mfbim))
      mdimtbx=mdimtbx+mmb(mfbim)
!
      enddo
!
      return
      end
end module
