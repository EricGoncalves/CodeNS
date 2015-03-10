module mod_clprd
implicit none
contains
      subroutine clprd( &
                 mfb,pres, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,l, &
                 pression,temp,cson)
!
!***********************************************************************
!
!_DA  DATE_C : avril 2002 - Eric GONCALVES / SINUMEF
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement de la condition aux limites : p=pres.
!_A    Utilisation des quatre relations de compatibilite discretisees
!_A    associee aux quatre caracteristiques sortantes
!_A      Dut=0
!_A      Dut=0
!_A      Dp - c2*Drho = 0
!_A      Dp + rho*c*Dun = 0
!_A    Normales interieures.
!
!     INP
!_I    ip11       : arg int              ; dim, nbr max de cellules de tous les
!_I                                        dom (pts fictifs inclus)
!_I    ip41       : arg int              ; dim, nbr max de pts de ttes les front
!_I    ip42       : arg int              ; dim, nbr max de pts de ttes les front
!_I                                        a normales stockees
!_I    ip60       : arg int              ; dim, nbr max d'equations
!_I    mfb        : arg int              ; numero de frontiere
!_I    pres       : arg real(ip41      ) ; pression statique a imposer
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    u          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    gam1       : com real             ; rap chal spec -1
!_I    gam4       : com real             ; 1/(rap chal spec -1)
!_I    gam5       : com real             ; rap chal spec/(rap chal spec -1)
!
!     COM
!_C    Notation 0      : valeurs a l' instant n.
!_C    Notation s      : valeurs issues du schema.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use proprieteflu
implicit none
integer :: mfb
double precision :: pres
integer :: ncbd
double precision :: v
integer :: mmb
integer :: mpb
integer :: mpn
integer :: l
double precision :: pression
double precision :: temp
double precision :: cson
double precision :: dqn
integer :: m
integer :: mb
integer :: mn
integer :: mt
integer :: n0c
integer :: n0n
integer :: nl
double precision :: ps
double precision :: qx
double precision :: qxs
double precision :: qy
double precision :: qys
double precision :: qz
double precision :: qzs
double precision :: rho
double precision :: roc0
!
!-----------------------------------------------------------------------
!
      real nxn,nyn,nzn
      dimension pres(ip40)
      dimension v(ip11,ip60)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
      dimension mmb(mtt),mpb(mtt)
      dimension mpn(mtt)
      dimension pression(ip11),temp(ip11),cson(ip11)
!
      n0n=npn(l)
      n0c=npc(l)
      mt=mmb(mfb)
!
      do m=1,mt
       mb=mpb(mfb)+m
       mn=mpn(mfb)+m
       nl=ncbd(mb)
!
       rho=v(nl,1)
       qxs=v(nl,2)/v(nl,1)
       qys=v(nl,3)/v(nl,1)
       qzs=v(nl,4)/v(nl,1)
       ps=gam1*(v(nl,5)-0.5*rho*(qxs**2+qys**2+qzs**2)-pinfl)
       roc0=sqrt(gam*ps*rho)
!
       dqn=(pres(m)-ps)/roc0
       qx=qxs+dqn*nxn(mn)   !attention aux signes des normales
       qy=qys+dqn*nyn(mn)
       qz=qzs+dqn*nzn(mn)
!
       v(nl,1)=v(nl,1)+rho*rho*dqn/roc0
       v(nl,2)=v(nl,1)*qx
       v(nl,3)=v(nl,1)*qy
       v(nl,4)=v(nl,1)*qz
       v(nl,5)=pres(m)/gam1+pinfl+0.5*v(nl,1)*(qx**2+qy**2+qz**2)
!
       pression(nl)=pres(m)
       temp(nl)=gam*pression(nl)/v(nl,1)
       cson(nl)=sqrt(temp(nl))
      enddo
!
      return
      end
end module
