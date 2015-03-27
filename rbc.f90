module mod_rbc
  use mod_rfvr
  use mod_rfvc
  implicit none
contains
  subroutine rbc( &
       ncin,nxn,nyn,nzn,ncbd, &
       sn,vol,v,mut, &
       bceqt, &
       rpi,rti,d0x,d0y,d0z,qtx,qty,qtz,x,y,z,omg, &
       pres,tp,rod,roud,rovd,rowd,roed, &
       icyc, &
       mnr,xnr,ynr,znr,mnc, &
       tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10,tm11, &
       tm12,tm13,pression,ztemp,cson)
!
!***********************************************************************
!
!     ACT
!_A    Remplissage des valeurs aux centres d'une rangee de mailles fictives
!_A    tout autour des domaines, par les valeurs aux centres des facettes
!_A    frontieres (conditions aux limites) et pour les domaines structures
!_A    les valeurs sur les aretes des maillages (extrapolations).
!
!     INP
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    sn         : arg real(ip31*ndir ) ; vecteur normal a une facette et de
!_I                                        norme egale a la surface de celle-ci
!_I    vol        : arg real(ip11      ) ; volume d'une cellule
!_I    u          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    rpi        : arg real(ip40      ) ; pres d'arret/pres d'arret etat de
!_I                                        ref utilisateur a imposer
!_I    rti        : arg real(ip40      ) ; temp d'arret/temp d'arret etat de
!_I                                        ref utilisateur a imposer
!_I    d0x        : arg real(ip40      ) ; composante en x d'une
!_I                                        direction de l'ecoulement
!_I    d0y        : arg real(ip40      ) ; composante en y d'une
!_I                                        direction de l'ecoulement
!_I    d0z        : arg real(ip40      ) ; composante en z d'une
!_I                                        direction de l'ecoulement
!_I    qtx        : arg real(ip40      ) ; composante en x d'une
!_I                                        vitesse tangentielle
!_I    qty        : arg real(ip40      ) ; composante en y d'une
!_I                                        vitesse tangentielle
!_I    qtz        : arg real(ip40      ) ; composante en z d'une
!_I                                        vitesse tangentielle
!_I    pres       : arg real(ip40      ) ; pression statique a imposer
!_I    tp         : arg real(ip40      ) ; temperature paroi a imposer
!_I    x          : arg real(ip21      ) ; coordonnee sur l'axe x
!_I    y          : arg real(ip21      ) ; coordonnee sur l'axe y
!_I    z          : arg real(ip21      ) ; coordonnee sur l'axe z
!_I    pres       : arg real(ip41      ) ; pression statique
!_I    it         : arg int              ; cycle courant du calcul
!_I    mnr        : arg int (ip44      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule recouvrante
!_I    xnr        : arg real(ip44      ) ; coefficient d'interpolation en x
!_I                                        dans une cellule recouvrante
!_I    ynr        : arg real(ip44      ) ; coefficient d'interpolation en y
!_I                                        dans une cellule recouvrante
!_I    znr        : arg real(ip44      ) ; coefficient d'interpolation en z
!_I                                        dans une cellule recouvrante
!_I    mnc        : arg int (ip43      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule coincidente
!_I    lzx        : com int              ; nbr total de domaines
!_I    mtnx       : com int              ; nbr total de frontieres
!_I                                        a normales stockes
!_I    mtcx       : com int              ; nbr total de frontieres coincidentes
!_I    mtrx       : com int              ; nbr total de frontieres recouvertes
!_I    mtax       : com int              ; nbr total de frontieres autres
!_I    td         : com char(lz        ) ; type de domaine (struct./non struct.)
!_I    nfbn       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une front a normales stockees
!_I    nfbc       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere coincidente
!_I    nfbr       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere recouverte
!_I    nfba       : com int (mtb       ) ; numero dans numerotation interne
!_I                                        d'une frontiere autre
!
!     I/O
!_/    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
!
!     LOC
!_L    nbd        : com int              ; nombre de frontieres a traiter
!_L    lbd        : com int (mtt       ) ; numero de front a traiter
!
!-----parameters figes--------------------------------------------------
!
    use para_var
    use para_fige
    use boundary
    use maillage
    use chainecarac
    use schemanum
    use mod_rbord
    use mod_rbse
    use mod_rfspstc
    use mod_rfve
    implicit none
    integer          ::       icyc,         l,        mf,       mfc,       mfr
    integer          ::  mnc(ip43), mnr(ip44),ncbd(ip41),ncin(ip41)
    double precision :: bceqt(ip41,neqt),      cson(ip11),       d0x(ip40),       d0y(ip40),       d0z(ip40)
    double precision ::        mut(ip12),       nxn(ip42),       nyn(ip42),       nzn(ip42),             omg
    double precision ::       pres(ip40),  pression(ip11),       qtx(ip40),       qty(ip40),       qtz(ip40)
    double precision ::        rod(ip40),      roed(ip40),      roud(ip40),      rovd(ip40),      rowd(ip40)
    double precision ::        rpi(ip40),       rti(ip40),   sn(ip31*ndir),       tm1(ip40),      tm10(ip40)
    double precision ::       tm11(ip40),      tm12(ip40),      tm13(ip40),       tm2(ip40),       tm3(ip40)
    double precision ::        tm4(ip40),       tm5(ip40),       tm6(ip40),       tm7(ip40),       tm8(ip40)
    double precision ::        tm9(ip40),        tp(ip40),    v(ip11,ip60),       vol(ip11),         x(ip21)
    double precision ::        xnr(ip44),         y(ip21),       ynr(ip44),         z(ip21),       znr(ip44)
    double precision ::      ztemp(ip11)
!
!-----------------------------------------------------------------------
!
!
!-----remplissage des valeurs de v aux points fictifs pour
!     les CL couplees avec le champs
!
    do mf=1,mtnx
       lbd(mf)=nfbn(mf)
    enddo
    nbd=mtnx
    call rfve( &
         v,pression,ztemp,cson, &
         ncbd,ncin)
!
    do mf=1,mtax
       lbd(mf)=nfba(mf)
    enddo
    nbd=mtax
    call rfve( &
         v,pression,ztemp,cson, &
         ncbd,ncin)
!
    do mfc=1,mtcx
       lbd(mfc)=nfbc(mfc)
    enddo
    nbd=mtcx
    call rfvc( &
         v,ncbd,mnc, &
         pression,ztemp,cson)
!
    do mfr=1,mtrx
       lbd(mfr)=nfbr(mfr)
    enddo
    nbd=mtrx
    call rfvr(v,ncbd,mnr,xnr,ynr,znr)
!
    call rbord( &
         1,1, &
         mnc,ncin,mnr,xnr,ynr,znr,nxn,nyn,nzn,ncbd, &
         sn,vol,v,vol, &
         bceqt, &
         rpi,rti,d0x,d0y,d0z,x,y,z, &
         pres,tp,rod,roud,rovd,rowd,roed, &
         icyc, &
         tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10,tm11, &
         tm12,tm13, &
         pression,ztemp,cson)
!
    if (equat(1:2).eq.'ns') then
       do mf=1,mtbx
          lbd(mf)=mf
       enddo
       nbd=mtbx
       call rbse(mut,ncin,ncbd)
    endif
    if (equat(6:7).eq.'ke') then
       call rbse(v(1,6),ncin,ncbd)
       call rbse(v(1,7),ncin,ncbd)
    endif
!
    do l=1,lzx
       call rfspstc(l,v(1,1))
       call rfspstc(l,v(1,2))
       call rfspstc(l,v(1,3))
       call rfspstc(l,v(1,4))
       call rfspstc(l,v(1,5))
       if(equat(1:2).eq.'ns') then
          call rfspstc(l,mut)
       endif
       if(equat(6:7).eq.'ke') then
          call rfspstc(l,v(1,6))
          call rfspstc(l,v(1,7))
       endif
    enddo
!
    return
  end subroutine rbc
end module mod_rbc
