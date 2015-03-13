module mod_cldebit
implicit none
contains
      subroutine cldebit( &
                 mfb,pres, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,l, &
                 pression,temp,cson)
!
!***********************************************************************
!
!_DA  DATE_C : janv 2006 - Eric Goncalves / LEGI
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement de la condition aux limites : debit impose.
!_A    Utilisation des quatre relations de compatibilite discretisees
!_A    associee aux quatre caracteristiques sortantes
!_A      Dut=0
!_A      Dut=0
!_A      Dp - c2*Drho = 0
!_A      Dp + rho*c*Dun = 0
!_A    Normales interieures.
!_A    ATTENTION: le debit local est stocke dans le tableau 'pres'.
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
double precision :: a2
double precision :: ca
double precision :: cb
double precision :: cc
double precision :: dp
double precision :: dqn
double precision :: drho
integer :: m
integer :: mb
integer :: mn
integer :: mt
integer :: nl
double precision :: ps
double precision :: qns
double precision :: qx
double precision :: qxs
double precision :: qy
double precision :: qys
double precision :: qz
double precision :: qzs
double precision :: rhos
double precision :: roc0
!
!-----------------------------------------------------------------------
!
      double precision nxn,nyn,nzn
!
      dimension pres(ip40)
      dimension v(ip11,ip60)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
      dimension mmb(mtt),mpb(mtt),mpn(mtt)
      dimension pression(ip11),temp(ip11),cson(ip11)
!
      mt=mmb(mfb)
      do m=1,mt
       mb=mpb(mfb)+m
       mn=mpn(mfb)+m
       nl=ncbd(mb)
!      etat schema
       rhos=v(nl,1)
       qxs=v(nl,2)/rhos
       qys=v(nl,3)/rhos
       qzs=v(nl,4)/rhos
       qns=qxs*nxn(mn)+qys*nyn(mn)+qzs*nzn(mn)
       ps=gam1*(v(nl,5)-0.5*rhos*(qxs**2+qys**2+qzs**2)-pinfl)
       roc0=sqrt(gam*ps*rhos)
       a2=gam*ps/rhos
!
       ca=a2
       cb=-ca*rhos+roc0*qns
       cc=roc0*pres(m)
       drho=0.5*(-cb+sqrt(abs(cb**2-4.*ca*cc)))/ca-rhos
       dp=a2*drho
       dqn=dp/roc0
       qx=qxs+dqn*nxn(mn)   !attention aux signes des normales
       qy=qys+dqn*nyn(mn)
       qz=qzs+dqn*nzn(mn)
!
       v(nl,1)=v(nl,1)+drho
       v(nl,2)=v(nl,1)*qx
       v(nl,3)=v(nl,1)*qy
       v(nl,4)=v(nl,1)*qz
       pression(nl)=ps+dp
       v(nl,5)=pression(nl)/gam1+pinfl+0.5*v(nl,1)*(qx**2+qy**2+qz**2)
       temp(nl)=gam*pression(nl)/v(nl,1)
       cson(nl)=sqrt(temp(nl))
      enddo
!
      return
      end subroutine
end module
