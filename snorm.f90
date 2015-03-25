module mod_snorm
  implicit none
contains
  subroutine snorm( &
       l,x,y,z, &
       sn,lgsnlt)
!
!***********************************************************************
!
!     ACT
!_A    Determination des surfaces orientees de toutes les facettes des
!_A    cellules d'un domaine structure.
!
!     INP
!_I    l          : arg int              ; numero de domaine
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    lgsnlt     : arg int              ; nombre de noeuds du dom (dont fic.)
!_I    equat      : com char             ; type d'equations modelisant l'ecoulement
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    ii2        : com int (lt        ) ; indice max en i reel
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jj2        : com int (lt        ) ; indice max en j reel
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!
!     OUT
!_O    sn         : arg real(lgsnlt,
!_O                          nind,ndir ) ; vecteur normal a une facette et de
!_O                                        norme egale a la surface de celle-ci
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use chainecarac
    use mod_norm
    implicit none
  integer          ::     i1,    i2,  imax,  imin,    j1
  integer          ::     j2,  jmax,  jmin,    k1,    k2
  integer          ::   kmax,  kmin,     l,lgsnlt
  double precision :: sn(lgsnlt,nind,ndir),             x(ip21),             y(ip21),             z(ip21)
!
!-----------------------------------------------------------------------
!
    character(len=7 ) :: eqt
!
    i1=ii1(l)
    i2=ii2(l)
    j1=jj1(l)
    j2=jj2(l)
    k1=kk1(l)
    k2=kk2(l)
!
    imin=i1-1
    imax=i2+1
    jmin=j1-1
    jmax=j2+1
    kmin=k1-1
    kmax=k2+1
    eqt=equat
!
    call norm( &
         l,x,y,z, &
         eqt,lgsnlt, &
         imin,imax,jmin,jmax,kmin,kmax, &
         sn(1,1,1),sn(1,1,2),sn(1,1,3), &
         sn(1,2,1),sn(1,2,2),sn(1,2,3), &
         sn(1,3,1),sn(1,3,2),sn(1,3,3))
!
    return
  end subroutine snorm
end module mod_snorm
