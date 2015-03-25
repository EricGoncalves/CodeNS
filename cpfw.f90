module mod_cpfw
  implicit none
contains
  subroutine cpfw( &
       ncyc, &
       x,y,z,r,exs1,exs2,nxn,nyn,nzn, &
       sn, &
       vol, &
       tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
       mu,mut,dist,cfke, &
       mnpar,fgam,utau, &
       v,dt, &
       ptdual,vdual,vdual1,vdual2, &
       tnte1,tnte2,tnte3,tnte4vv, &
       toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
       tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10, &
       tm11,tm12,tm13, &
       ncin, &
       mnc, &
       ncbd,mnr,xnr,ynr,znr, &
       bceqt, &
       rpi,rti,d0x,d0y,d0z,qtx,qty,qtz,pres,tp, &
       rod,roud,rovd,rowd,roed, &
       pression,ztemp,cson, &
       cvi,cvj,cvk, &
       cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
!***********************************************************************
!
!     ACT
!_A    Resolution des equations d'Euler ou de Navier-Stokes moyennees
!_A    completees par un modele de turbulence par une methode
!_A    (pseudo-)instationnaire pendant ncycle(l) cycles
!_A    sur chaque niveau de grille l specifie.
!
!     INP
!_I    ncycle     : arg int              ; nbr tot de cycles de l'execution courante
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    exs1       : arg real             ; premier coef d'interpolation a l'
!_I                                        l'ordre 0 du couple de coef exs1,exs2
!_I    exs2       : arg real             ; deuxieme coef d'interpolation a l'
!_I                                        l'ordre 0 du couple de coef exs1,exs2
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    mnc        : arg int (ip43      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule coincidente
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mnr        : arg int (ip44      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule recouvrante
!_I    xnr        : arg real(ip44      ) ; coefficient d'interpolation en x
!_I                                        dans une cellule recouvrante
!_I    ynr        : arg real(ip44      ) ; coefficient d'interpolation en y
!_I                                        dans une cellule recouvrante
!_I    znr        : arg real(ip44      ) ; coefficient d'interpolation en z
!_I                                        dans une cellule recouvrante
!_I    rpi        : arg real(ip40      ) ; pres d'arret/pres d'arret etat de
!_I                                        ref utilisateur a imposer
!_I    rti        : arg real(ip40      ) ; temp d'arret/temp d'arret etat de
!_I                                        ref utilisateur a imposer
!_I    d0x        : arg real(ip40      ) ; composante en x d'une direction
!_I                                        de l'ecoulement a imposer
!_I    d0y        : arg real(ip40      ) ; composante en y d'une direction
!_I                                        de l'ecoulement a imposer
!_I    d0z        : arg real(ip40      ) ; composante en z d'une direction
!_I                                        de l'ecoulement a imposer
!_I    qtx        : arg real(ip40      ) ; composante en x d'une
!_I                                        vitesse tangentielle a imposer
!_I    qty        : arg real(ip40      ) ; composante en y d'une
!_I                                        vitesse tangentielle a imposer
!_I    qtz        : arg real(ip40      ) ; composante en z d'une
!_I                                        vitesse tangentielle a imposer
!_I    pres       : arg real(ip40      ) ; pression statique a imposer
!_I    tp         : arg real(ip40      ) ; temperature paroi a imposer
!_I    titrt1     : com char             ; titre du calcul
!_I    equat      : com char             ; type d'equations modelisant l'ecoulement
!_I    out        : com int              ; unite logiq, moyennes des residus
!_I    imp        : com int              ; unite logiq, sorties de controle
!_I    lzx        : com int              ; nbr total de domaines
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    mtnx       : com int              ; nbr total de frontieres a normales stockes
!_I    mtcx       : com int              ; nbr total de frontieres coincidentes
!_I    mtrx       : com int              ; nbr total de frontieres recouvertes
!_I    mtax       : com int              ; nbr total de frontieres autres
!_I    nnn        : com int (lt        ) ; nombre de noeuds du dom (dont fic.)
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
!_I    nnc        : com int (lt        ) ; nombre de cellules du dom (dont fic.)
!_I    npfb       : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes facettes
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    ii1        : com int (lt        ) ; indice min en i reel
!_I    ii2        : com int (lt        ) ; indice max en i reel
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jj1        : com int (lt        ) ; indice min en j reel
!_I    jj2        : com int (lt        ) ; indice max en j reel
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    kd1        : com int (lt        ) ; indice min en k fictif
!_I    kk1        : com int (lt        ) ; indice min en k reel
!_I    kk2        : com int (lt        ) ; indice max en k reel
!_I    kd2        : com int (lt        ) ; indice max en k fictif
!_I    mmb        : com int (mtt       ) ; nombre de pts d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    iminb      : com int (mtt       ) ; indice min en i d'une frontiere
!_I    imaxb      : com int (mtt       ) ; indice max en i d'une frontiere
!_I    jminb      : com int (mtt       ) ; indice min en j d'une frontiere
!_I    jmaxb      : com int (mtt       ) ; indice max en j d'une frontiere
!_I    kminb      : com int (mtt       ) ; indice min en k d'une frontiere
!_I    kmaxb      : com int (mtt       ) ; indice max en k d'une frontiere
!_I    nfbr       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere recouverte
!_I    nfbc       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere coincidente
!_I    nfbn       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une front a normales stockees
!_I    kimp       : com int              ; niveau de sortie sur unite logi imp
!_I    klomg      : com int              ; cle pour rotation du repere relatif
!_I    omg        : com real             ; vitesse rotation du repere relatif
!_I    icytur0    : com int              ; nbr de cycl en deb de calcul au cours
!_I                                        desquelles mut n'est pas mis a jour
!_I    ncyturb    : com int              ; freq en it de mise a jour de mut
!_I    kdtl       : com int              ; cle d'utilisation pas de temps local
!_I    icychr0    : com int              ; nbr de cycl en deb de calcul au cours
!_I                                        desquelles le pas de temps mis a jour
!_I                                        a chaque it
!_I    ncychro    : com int              ; freq en it de mise a jour du
!_I                                        pas de temps
!_I    dt1min     : com real             ; pas de temps constant desire
!_I    eta        : com real(lt        ) ; nombre de CFL
!_I    ncyresi    : com int              ; freq en it de calcul des residus
!_I    ncysave    : com int              ; freq en it de sauvegarde des var aero
!_I    kmf        : com int (lt        ) ; cle phase implicite
!
!     OUT
!_O    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_O                                        norme egale a la surface de celle-ci
!_O    vol        : arg real(ip11      ) ; volume d'une cellule
!_O    mu         : arg real(ip12      ) ; viscosite moleculaire
!_O    dist       : arg real(ip12      ) ; distance a la paroi
!_O    nfba       : com int (mtb       ) ; numero dans numerotation interne
!_O                                        d'une frontiere autre
!
!     I/O
!_/    numt       : arg int              ; cycle courant du calcul
!_/    ncyc       : arg int              ; cycle courant de l'execution
!_/    mut        : arg real(ip12      ) ; viscosite turbulente
!_/    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!
!     LOC
!_L    r          : arg real(ip11      ) ; distance a l'axe
!_L    u          : arg real(ip11,ip60 ) ; variables a l'instant n
!_L    dt         : arg real(ip11      ) ; pas de temps
!_L    d          : arg real(ip11,ip60 ) ; terme de dissipation artificielle
!_L    toxx       : arg real(ip12      ) ; composante en xx du tenseur des
!_L                                        contraintes visqueuses
!_L    toxy       : arg real(ip12      ) ; composante en xy du tenseur des
!_L                                        contraintes visqueuses
!_L    toxz       : arg real(ip12      ) ; composante en xz du tenseur des
!_L                                        contraintes visqueuses
!_L    toyy       : arg real(ip12      ) ; composante en yy du tenseur des
!_L                                        contraintes visqueuses
!_L    toyz       : arg real(ip12      ) ; composante en yz du tenseur des
!_L                                        contraintes visqueuses
!_L    tozz       : arg real(ip12      ) ; composante en zz du tenseur des
!_L                                        contraintes visqueuses
!_L    qcx        : arg real(ip12      ) ; composante en x du flux de chaleur
!_L    qcy        : arg real(ip12      ) ; composante en y du flux de chaleur
!_L    qcz        : arg real(ip12      ) ; composante en z du flux de chaleur
!_L    nbd        : com int              ; nombre de frontieres a traiter
!_L    lbd        : com int (mtt       ) ; numero de front a traiter
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use sortiefichier
    use boundary
    use chainecarac
    use schemanum
    use maillage
    use definition
    use modeleturb
    use mod_smg_num
    use mod_metric
    use mod_metric3
    implicit none
  integer          ::         img,          l,        lev,         lm,        mfb
  integer          ::         mfc,        mfn,        mfr,  mnc(ip43),mnpar(ip12)
  integer          ::   mnr(ip44), ncbd(ip41), ncin(ip41),       ncyc,        ngx
  double precision ::   bceqt(ip41,neqt),        cfke(ip13),       cmui1(ip21),       cmui2(ip21),       cmuj1(ip21)
  double precision ::        cmuj2(ip21),       cmuk1(ip21),       cmuk2(ip21),        cson(ip11),         cvi(ip21)
  double precision ::          cvj(ip21),         cvk(ip21),         d0x(ip40),         d0y(ip40),         d0z(ip40)
  double precision ::         dist(ip12),          dt(ip11),              exs1,              exs2,        fgam(ip42)
  double precision ::           mu(ip12),         mut(ip12),         nxn(ip42),         nyn(ip42),         nzn(ip42)
  double precision ::         pres(ip40),    pression(ip11), ptdual(ip11,ip60),         qcx(ip12),         qcy(ip12)
  double precision ::          qcz(ip12),         qtx(ip40),         qty(ip40),         qtz(ip40),           r(ip11)
  double precision ::          rod(ip40),        roed(ip40),        roud(ip40),        rovd(ip40),        rowd(ip40)
  double precision ::          rpi(ip40),         rti(ip40),     sn(ip31*ndir),         tm1(ip40),        tm10(ip40)
  double precision ::         tm11(ip40),        tm12(ip40),        tm13(ip40),         tm2(ip40),         tm3(ip40)
  double precision ::          tm4(ip40),         tm5(ip40),         tm6(ip40),         tm7(ip40),         tm8(ip40)
  double precision ::          tm9(ip40),         tn1(ip00),        tn10(ip00),         tn2(ip00),         tn3(ip00)
  double precision ::          tn4(ip00),         tn5(ip00),         tn6(ip00),         tn7(ip00),         tn8(ip00)
  double precision ::          tn9(ip00),  tnte1(ip11,ip60),  tnte2(ip11,ip60),  tnte3(ip11,ip60),tnte4vv(ip11,ip60)
  double precision ::         toxx(ip12),        toxy(ip12),        toxz(ip12),        toyy(ip12),        toyz(ip12)
  double precision ::         tozz(ip12),          tp(ip40),        utau(ip42),      v(ip11,ip60),  vdual(ip11,ip60)
  double precision ::  vdual1(ip11,ip60), vdual2(ip11,ip60),         vol(ip11),           x(ip21),         xnr(ip44)
  double precision ::            y(ip21),         ynr(ip44),           z(ip21),         znr(ip44),       ztemp(ip11)
!
!-----------------------------------------------------------------------
!
    character(len=1316) :: form
!
!
!-----calcul de la metrique-------------------------------------------
!
    if(ncyc.eq.0) then
!
       call metric( &
            x,y,z,r,exs1,exs2, &
            sn,vol, &
            ncbd,mnc, &
            mnr,xnr,ynr,znr, &
            tn1,tn2,tn3)

       if((ischema.eq.2).or.(ischema.eq.3).or.(ischema.eq.4).or.(ischema.eq.6).or. &
            (ischema.eq.8).or.(ischema.eq.11).or.(ischema.eq.13).or.(ischema.eq.15)) then
          do l=1,lzx
             do  img=1,lgx
                lm=l+(img-1)*lz
                call metric3( &
                     lm,x,y,z, &
                     cvi,cvj,cvk, &
                     cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
             enddo
          enddo
       endif
!        do l=1,lzx
!         do  img=1,lgx
!          lm=l+(img-1)*lz
!          call metric2( &
!                 lm,x,y,z, &
!                 cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!         enddo
!        enddo
!
!-----reperage des frontieres autres------------------------------------
!
       do mfb=1,mtbx
          lbd(mfb)=1
       enddo
!
       do mfn=1,mtnx
          mfb=nfbn(mfn)
          lbd(mfb)=0
       enddo
!
       do mfc=1,mtcx
          mfb=nfbc(mfc)
          lbd(mfb)=0
       enddo
!
       do mfr=1,mtrx
          mfb=nfbr(mfr)
          lbd(mfb)=0
       enddo
!
       do mfb=1,mtbx
          if(lbd(mfb).eq.1) then
             mtax=mtax+1
             nfba(mtax)=mfb
          endif
       enddo
!
    endif  !fin test sur ncyc
!
    ngx=1
    if(kfmg.ne.0.and.kfmg.ne.3) ngx=lgx
!
!-----boucle iterative sur les niveaux de grille-----------------
!
    do lev=1,ngx
       img = ngx - lev + 1
       mgl = img
       write(imp,*)' '
       write(imp,*)' ITERATION POUR LA GRILLE =',img
       write(imp,*)' Nombre de cycles requis =',ncycle(img)
!
       if(ncycle(img).ne.0)then
!
          call smg_num( &
               img,ncycle(img), &
               numt, &
               x,y,z,r,nxn,nyn,nzn, &
               sn, &
               vol, &
               tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10, &
               mu,mut,dist,cfke, &
               mnpar,fgam,utau, &
               v,dt, &
               ptdual,vdual,vdual1,vdual2, &
               tnte1,tnte2,tnte3,tnte4vv, &
               toxx,toxy,toxz,toyy,toyz,tozz,qcx,qcy,qcz, &
               tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10, &
               tm11,tm12,tm13, &
               ncin, &
               mnc, &
               ncbd,mnr,xnr,ynr,znr, &
               bceqt, &
               rpi,rti,d0x,d0y,d0z,qtx,qty,qtz,pres,tp, &
               rod,roud,rovd,rowd,roed, &
               pression,ztemp,cson, &
               cvi,cvj,cvk, &
               cmui1,cmui2,cmuj1,cmuj2,cmuk1,cmuk2)
!
       endif
!
!     mise a jour de variables a l'interface du raccord
!     -------------------------------------------------
!
       do mfc=1,mtcx
          lbd(mfc)=nfbc(mfc)+(img-1)*mtb
       enddo
       nbd=mtcx
       call rfvc( & 
            v,ncbd,mnc, &
            pression,ztemp,cson)
!            call rfvc(v,ncbd,mnc)
!
       do mfr=1,mtrx
          lbd(mfr)=nfbr(mfr)+(img-1)*mtb
       enddo
       nbd=mtrx
       call rfvr( &
            v, &
            ncbd,mnr,xnr,ynr,znr)
!
!     cell-to-node projection of the variables at level "img"
!
       call smg_cn( &
            img, &
            vol,tnte2(1,1), &
            v,tnte4vv)
!
!     node-to-cell transfer of the variables from level "img" -->"img-1"
!
       if(lev.lt.ngx) then
!
          call smg_cf( &
               img,img-1, &
               vol, &
               tnte4vv,v)
!
       endif
    enddo !fin des iterations
!
    return
  end subroutine cpfw
end module mod_cpfw
