module mod_zvisqc
  implicit none
contains
  subroutine zvisqc( &
       img, &
       s,mu,ro, &
       x,y,z,mut,dist,mnpar,fgam, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       icyc,mcyturb, &
       ncbd,ncin,mnc, &
       sn,vol, &
       dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2, &
       ztemp)
!
!***********************************************************************
!
!     ACT
!_A    Determination du tenseur des contraintes visqueuses et
!_A    du flux de chaleur.
!
!_I    s          : arg real(ip11,ip60 ) ; variables de calcul
!_I    mu         : arg real(ip12      ) ; viscosite moleculaire
!_I    ro         : arg real(ip11      ) ; masse volumique
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    r          : arg real(ip11      ) ; distance a l'axe
!_I    mut        : arg real(ip12      ) ; viscosite turbulente
!_I    dist       : arg real(ip12      ) ; distance a la paroi
!_I    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_I                                        contraintes visqueuses
!_I    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_I                                        contraintes visqueuses
!_I    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_I                                        contraintes visqueuses
!_I    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_I                                        contraintes visqueuses
!_I    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_I                                        contraintes visqueuses
!_I    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_I                                        contraintes visqueuses
!_I    qcx        : arg real(ip12      ) ; composante en x du flux de chaleur
!_I    qcy        : arg real(ip12      ) ; composante en y du flux de chaleur
!_I    qcz        : arg real(ip12      ) ; composante en z du flux de chaleur
!_I    it         : arg int              ; cycle courant du calcul
!_I    mcyturb    : arg int              ; reste de division de nit par ncyturb
!_I    icytur0    : arg int              ; nbr de cycl en deb de calcul au cours
!_I                                        desquelles mut n'est pas mis a jour
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    mnc        : arg int (ip43      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule coincidente
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    equat      : com char             ; type d'equations modelisant l'ecoule-
!_I                                        ment
!_I    nnn        : com int (lt        ) ; nombre de noeuds du dom (dont fic.)
!_I    npfb       : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes facettes
!_I    nfbn       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une front a normales stockees
!_I    nfbc       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere coincidente
!_I    nfbr       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere recouverte
!_I    nfba       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere autre
!_I    lzx        : com int              ; nbr total de domaines
!_I    mtnx       : com int              ; nbr total de frontieres
!_I                                        a normales stockes
!_I    mtcx       : com int              ; nbr total de frontieres coincidentes
!_I    mtrx       : com int              ; nbr total de frontieres recouvertes
!_I    mtax       : com int              ; nbr total de frontieres autres
!
!     LOC
!_L    vxtn1      : arg real(ip00      ) ; comp en x de la vit ou tab de travail
!_L    vy         : arg real(ip00      ) ; composante en y de la vitesse
!_L    vz         : arg real(ip00      ) ; composante en z de la vitesse
!_L    ztemp      : arg real(ip00      ) ; temperature
!_L    dvxx       : arg real(ip00      ) ; composante xx du gradient de vitesse
!_L    dvxy       : arg real(ip00      ) ; composante xy du gradient de vitesse
!_L    dvxz       : arg real(ip00      ) ; composante xz du gradient de vitesse
!_L    dvyx       : arg real(ip00      ) ; composante yx du gradient de vitesse
!_L    dvyy       : arg real(ip00      ) ; composante yy du gradient de vitesse
!_L    dvyz       : arg real(ip00      ) ; composante yz du gradient de vitesse
!_L    dvzx       : arg real(ip00      ) ; composante zx du gradient de vitesse
!_L    dvzy       : arg real(ip00      ) ; composante zy du gradient de vitesse
!_L    dvzz       : arg real(ip00      ) ; composante zz du gradient de vitesse
!_L    dtx        : arg real(ip00      ) ; composante x du gradient de temp
!_L    dty        : arg real(ip00      ) ; composante y du gradient de temp
!_L    dtz        : arg real(ip00      ) ; composante z du gradient de temp
!_L    nbd        : com int              ; nombre de frontieres a traiter
!_L    lbd        : com int (mtt       ) ; numero de front a traiter
!
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    use maillage
    use chainecarac
    use schemanum
    use modeleturb
    use mod_smg_fcs
    use mod_zfluto
    use mod_atctranske
    use mod_zgrad
    use mod_rbte
    use mod_zgrad2
    use mod_rbtc
    implicit none
    integer          ::        icyc,        img,          l,     lgsnlt,         lm
    integer          ::     mcyturb,         mf,        mfc,  mnc(ip43),mnpar(ip12)
    integer          ::  ncbd(ip41), ncin(ip41),       npsn
    double precision ::   cmui1(ip21),  cmui2(ip21),  cmuj1(ip21),  cmuj2(ip21),  cmuk1(ip21)
    double precision ::   cmuk2(ip21),   dist(ip12),   dvxx(ip00),   dvxy(ip00),   dvxz(ip00)
    double precision ::    dvyx(ip00),   dvyy(ip00),   dvyz(ip00),   dvzx(ip00),   dvzy(ip00)
    double precision ::    dvzz(ip00),   fgam(ip42),     mu(ip12),    mut(ip12),    qcx(ip12)
    double precision ::     qcy(ip12),    qcz(ip12),     ro(ip11), s(ip11,ip60),sn(ip31*ndir)
    double precision ::    toxx(ip12),   toxy(ip12),   toxz(ip12),   toyy(ip12),   toyz(ip12)
    double precision ::    tozz(ip12),    vol(ip11),      x(ip21),      y(ip21),      z(ip21)
    double precision ::   ztemp(ip11)
    double precision,allocatable :: dtx(:),dty(:),dtz(:)
!
!-----------------------------------------------------------------------
!
    ALLOCATE(dtx(ip00),dty(ip00),dtz(ip00))
!
    do l=1,lzx
       lm=l+(img-1)*lz
       npsn  =ndir*npfb(lm)+1
       lgsnlt=nnn(lm)
!
       if(ischema.eq.1) then
          call zgrad( &
               lm, &
               equat, &
               sn(npsn),lgsnlt, &
               vol, &
               s,ztemp, &
               dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
               dtx,dty,dtz)
       else
          call zgrad2( &
               lm, &
               equat, &
               sn(npsn),lgsnlt, &
               vol, &
               s,ztemp, &
               dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
               dtx,dty,dtz, &
               cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
       endif
!
       if(icyc.gt.icytur0) then
          if(mcyturb.eq.1) then
             if(img.eq.mgl) then!
                if(ktransi.eq.1) then  ! transition fixee
                   call atctranske(l,s,mu,mut,mnpar,fgam)
                endif
             else  ! fine to coarse transfer of mut
                call smg_fcs( &
                     img-1,img, &
                     vol,mut)
             endif
          endif
       endif
!
       call zfluto( &
            lm,mu,mut,toxx,toxy,toxz,toyy,toyz,tozz, &
            qcx,qcy,qcz, &
            dvxx,dvxy,dvxz,dvyx,dvyy,dvyz,dvzx,dvzy,dvzz, &
            dtx,dty,dtz)
!
    enddo
!
    do mf=1,mtnx
       lbd(mf)=nfbn(mf)+(img-1)*mtb
    enddo
    nbd=mtnx
    call rbte( &
         toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
         ncbd,ncin)
!
    do mf=1,mtax
       lbd(mf)=nfba(mf)+(img-1)*mtb
    enddo
    nbd=mtax
    call rbte( &
         toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
         ncbd,ncin)
!
!
    do mfc=1,mtcx
       lbd(mfc)=nfbc(mfc)+(img-1)*mtb
    enddo
    nbd=mtcx
    call rbtc( &
         toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
         ncbd,ncin,mnc)

    DEALLOCATE(dtx,dty,dtz)

    return
  end subroutine zvisqc
end module mod_zvisqc
