      subroutine clgli2( &
                 ncbd,ncin, &
                 mfb,mmb,mpb,mpn, &
                 nxn,nyn,nzn, &
                 v,pression,temp,cson)
!
!***********************************************************************
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement de la condition de glissement (qn=0.)
!_A    Normales interieures.
!
!     INP
!_I    mfb        : arg int              ; numero de frontiere
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    gam1       : com real             ; gam -1
!_I    gam3       : com real             ; 1/gam
!_I    gam4       : com real             ; 1/(gam-1)
!_I    gam5       : com real             ; gam/(gam-1)
!
!     COM
!_C    Notation 0      : valeurs a l' instant n.
!_C    Notation s      : valeurs issues du schema.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use proprieteflu
implicit none
integer :: ncbd
integer :: ncin
integer :: mfb
integer :: mmb
integer :: mpb
integer :: mpn
double precision :: v
double precision :: pression
double precision :: temp
double precision :: cson
integer :: m
integer :: mb
integer :: mn
integer :: mt
integer :: ni
integer :: nl
double precision :: qn
double precision :: qtx
double precision :: qty
double precision :: qtz
double precision :: qx
double precision :: qy
double precision :: qz
double precision :: rho
!
!-----------------------------------------------------------------------
!
      real nxn,nyn,nzn
      dimension v(ip11,ip60)
      dimension ncbd(ip41),ncin(ip41)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
      dimension mmb(mtt),mpb(mtt),mpn(mtt)
      dimension pression(ip11),temp(ip11),cson(ip11)
!
      mt=mmb(mfb)
!
      do m=1,mt
       mb=mpb(mfb)+m
       mn=mpn(mfb)+m
       nl=ncbd(mb)
       ni=ncin(mb)
!
       rho=v(ni,1)
       qx=v(ni,2)/rho
       qy=v(ni,3)/rho
       qz=v(ni,4)/rho
       qn=qx*nxn(mn)+qy*nyn(mn)+qz*nzn(mn)
       qtx=qx-qn*nxn(mn)
       qty=qy-qn*nyn(mn)
       qtz=qz-qn*nzn(mn)
!
       v(nl,1)=rho
       v(nl,2)=v(nl,1)*qtx
       v(nl,3)=v(nl,1)*qty
       v(nl,4)=v(nl,1)*qtz
       pression(nl)=pression(ni)      !hypothese dPdn=0
       v(nl,5)=pression(nl)/gam1+pinfl+0.5*rho*(qtx**2+qty**2+qtz**2)
       temp(nl)=temp(ni)              !flux de chaleur nul: dTdn=0
       cson(nl)=cson(ni)
      enddo
!
      return
      end
