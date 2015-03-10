      subroutine chdgcv( &
                 l,disc, &
                 x,y,z, &
                 tn1,tn2,tn3, &
                 imin,imax,jmin,jmax,kmin,kmax,kdg)
!
!***********************************************************************
!
!     ACT
!_A    Changement eventuel de discretisation pour les coordonnees
!_A    a partir de valeurs connues aux noeuds du maillage.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    disc       : arg char             ; changement de discretisation (centre/noeud)
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    kdgv       : com int              ; unite logiq, coordonnees des noeuds
!_I    kdgc       : com int              ; unite logiq, coordonnees des centres
!_I    kdgcf      : com int              ; unite logiq, coordonnees des centres
!_I                                        et centres des facettes front
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    ii2        : com int (lt        ) ; indice max en i reel
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jj2        : com int (lt        ) ; indice max en j reel
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!
!     OUT
!_O    tn1        : arg real(ip00      ) ; coordonnee sur l'axe x
!_O    tn2        : arg real(ip00      ) ; coordonnee sur l'axe y
!_O    tn3        : arg real(ip00      ) ; coordonnee sur l'axe z
!_O    imin       : arg int              ; indice min en i
!_O    imax       : arg int              ; indice max en i
!_O    jmin       : arg int              ; indice min en j
!_O    jmax       : arg int              ; indice max en j
!_O    kmin       : arg int              ; indice min en k
!_O    kmax       : arg int              ; indice max en k
!_O    kdg        : arg int              ; unite logique, maillage
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
   use maillage
   use sortiefichier
implicit none
integer :: l
double precision :: x
double precision :: y
double precision :: z
double precision :: tn1
double precision :: tn2
double precision :: tn3
integer :: imin
integer :: imax
integer :: jmin
integer :: jmax
integer :: kmin
integer :: kmax
integer :: kdg
!
!-----------------------------------------------------------------------
!
      character(len=4 ) :: disc
      dimension x(ip21),y(ip21),z(ip21)
      dimension tn1(ip00),tn2(ip00),tn3(ip00)
!
      if(disc.eq.'cvcc') then
        call cvccg( &
                 l, &
                 x,y,z, &
                 tn1,tn2,tn3)
      imin=ii1(l)
      imax=ii2(l)-1
      jmin=jj1(l)
      jmax=jj2(l)-1
      kmin=kk1(l)
      kmax=kk2(l)-1
      kdg=kdgc
      elseif(disc.eq.'cvcf') then
            call cvccg( &
                 l, &
                 x,y,z, &
                 tn1,tn2,tn3)
      imin=ii1(l)-1
      imax=ii2(l)
      jmin=jj1(l)-1
      jmax=jj2(l)
      kmin=kk1(l)-1
      kmax=kk2(l)
      kdg=kdgcf
      elseif(disc.eq.'cvcv') then
            call cvcvg( &
                 l, &
                 x,y,z, &
                 tn1,tn2,tn3)
      imin=ii1(l)
      imax=ii2(l)
      jmin=jj1(l)
      jmax=jj2(l)
      kmin=kk1(l)
      kmax=kk2(l)
      kdg=kdgv
      endif
!
      return
      end
