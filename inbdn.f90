module mod_inbdn
  implicit none
contains
  subroutine inbdn( &
       mfbe,kibdn, &
       x,y,z, &
       sn, &
       ncbd,nxn,nyn,nzn, &
       tn1,tn2,tn3,tn4,tn5,tn6, &
       tn7,tn8,tn9)
!
!***********************************************************************
!
!     ACT
!_A    Initialisation des donnees necessaires aux frontieres a normales
!_A    stockees, principalement pour chaque point les trois composantes
!_A    de la normale a la facette frontiere.
!
!     INP
!_I    mfbe       : arg int              ; numero externe de frontiere
!_I    kibdn      : arg int              ; cle initialisat. normales a la front
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    indfl      : com char(mtb       ) ; type de plan de la frontiere
!_I    equat      : com char             ; type d'equations modelisant l'ecoule-ment
!_I    npfb       : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes facettes
!_I    nnfb       : com int (lt        ) ; nombre de facettes du dom (dont fic.)
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    nfei       : com int (mtb       ) ; numero de base interne d'une front
!_I                                        en fct du numero externe
!_I    iminb      : com int (mtt       ) ; indice min en i d'une frontiere
!_I    imaxb      : com int (mtt       ) ; indice max en i d'une frontiere
!_I    jminb      : com int (mtt       ) ; indice min en j d'une frontiere
!_I    jmaxb      : com int (mtt       ) ; indice max en j d'une frontiere
!_I    kminb      : com int (mtt       ) ; indice min en k d'une frontiere
!_I    kmaxb      : com int (mtt       ) ; indice max en k d'une frontiere
!_I    kfn        : com int              ; unite logiq, tableaux normales
!
!     OUT
!_O    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_O                                        normal a une facette frontiere
!_O    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_O                                        normal a une facette frontiere
!_O    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_O                                        normal a une facette frontiere
!_O    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_O                                        dans tab front a normales stockees
!_O    nfbn       : com int (mtb       ) ; numero dans numerotation interne
!_O                                        d'une front a normales stockees
!
!     I/O
!_/    mtnx       : com int              ; nbr total de frontieres
!_/                                        a normales stockes
!_/    mdimtnx    : com int              ; nbr de pts de ttes les front
!_/                                        a normales stockees
!
!     LOC
!_L    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_L                                        norme egale a la surface de celle-ci
!_L    tn1        : arg real(ip00      ) ; tableau de travail
!_L    tn2        : arg real(ip00      ) ; tableau de travail
!_L    tn3        : arg real(ip00      ) ; tableau de travail
!_L    tn4        : arg real(ip00      ) ; tableau de travail
!_L    tn5        : arg real(ip00      ) ; tableau de travail
!_L    tn6        : arg real(ip00      ) ; tableau de travail
!_L    tn7        : arg real(ip00      ) ; tableau de travail
!_L    tn8        : arg real(ip00      ) ; tableau de travail
!_L    tn9        : arg real(ip00      ) ; tableau de travail
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use maillage
    use boundary
    use chainecarac
    use sortiefichier
    use mod_initns
    implicit none
    integer          ::  imax,  img, imin, jmax, jmin
    integer          :: kibdn, kmax, kmin,    l,   lm
    integer          ::   m0n, mfbe, mfbi,mfbim,  mfn
    integer          ::    mt, ncbd
    double precision :: nxn,nyn,nzn, sn,tn1
    double precision :: tn2,tn3,tn4,tn5,tn6
    double precision :: tn7,tn8,tn9,  x,  y
    double precision ::   z
!
!-----------------------------------------------------------------------
!
    character(len=2 ) :: indfb
    character(len=7 ) :: eqt
    dimension x(ip21),y(ip21),z(ip21)
    dimension sn(ip31*ndir)
    dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
    dimension tn1(ip00),tn2(ip00),tn3(ip00), &
              tn4(ip00),tn5(ip00),tn6(ip00), &
              tn7(ip00),tn8(ip00),tn9(ip00)
!
    mfbi=nfei(mfbe)
!
    mtnx=mtnx+1
    mfn =mtnx
    nfbn(mfn)=mfbi
    l=ndlb(mfbi)
!
    do img=1,lgx
!
       mfbim=mfbi+(img-1)*mtb
       mt=mmb(mfbim)
       mpn(mfbim)=mdimtnx
       m0n=mpn(mfbim)
       mdimtnx=mdimtnx+mt
!
!     remplissage des tableaux  nxn , nyn , nzn
!
       if((kibdn.eq.1).or.(img.gt.1)) then
!
          lm=l+(img-1)*lz
!
          imin=iminb(mfbim)
          imax=imaxb(mfbim)
          jmin=jminb(mfbim)
          jmax=jmaxb(mfbim)
          kmin=kminb(mfbim)
          kmax=kmaxb(mfbim)
          indfb=indfl(mfbi)
          eqt=equat

          call initns( &
               mfbim,lm,indfb, &
               imin,imax,jmin,jmax,kmin,kmax, &
               eqt,x,y,z, &
               ncbd,nxn,nyn,nzn, &
               tn1,tn2,tn3,tn4,tn5,tn6, &
               tn7,tn8,tn9)
!
       elseif(kibdn.eq.0) then
!       call readfn( &
!                 kfn,nxn,nyn,nzn, &
!                 mt,m0n)
       endif

    enddo
!
    return
  end subroutine inbdn
end module mod_inbdn
