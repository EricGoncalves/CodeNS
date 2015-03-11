module mod_clprd_prcd
implicit none
contains
      subroutine clprd_prcd( &
                 mfb,pres, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,l, &
                 pression,temp,cson)
!
!***********************************************************************
!
!_DA   aout 2002 - AUTEUR: Eric GONCALVES / SINUMEF
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement de la condition aux limites : p=pres.
!_A    Utilisation de quatre relations de compatibilite discretisees
!_A    associee aux quatre caracteristiques sortantes
!_A      Dut=0    (vitesse tangentielle)
!_A      Dut=0
!_A      Dp - c2*Drho = 0
!_A      (lambda+ -un)*Dp + rho*beta2*c2*Dun = 0    (vitesse normale)
!_A    Normales interieures.
!_A    Preconditionnement basse vitesse de Turkel
!
!
!     INP
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
!
!     COM
!_C    Notation s      : valeurs issues du schema.
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use proprieteflu
      use definition
      use schemanum
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
double precision :: alp
double precision :: beta2
double precision :: dqn
integer :: m
integer :: mb
integer :: mn
integer :: mt
integer :: n0c
integer :: n0n
integer :: nl
double precision :: ps
double precision :: q2
double precision :: qns
double precision :: qx
double precision :: qxs
double precision :: qy
double precision :: qys
double precision :: qz
double precision :: qzs
double precision :: rho
double precision :: roc
!
!-----------------------------------------------------------------------
!
      real nxn,nyn,nzn,qinf
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
      qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
!
      do m=1,mt
       mb=mpb(mfb)+m
       mn=mpn(mfb)+m
       nl=ncbd(mb)
       rho=v(nl,1)
       qxs=v(nl,2)/rho
       qys=v(nl,3)/rho
       qzs=v(nl,4)/rho
       qns=qxs*nxn(mn)+qys*nyn(mn)+qzs*nzn(mn)
       q2=qxs**2+qys**2+qzs**2
       ps=gam1*(v(nl,5)-0.5*rho*q2-pinfl)
       roc=sqrt(gam*ps*rho)
       a2=gam*ps/rho
       beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
       alp=0.5*((1.+beta2)*abs(qns) + sqrt(((1.-beta2)*qns)**2+4.*beta2*a2))
!
       dqn=(alp+qns)*(pres(m)-ps)/(rho*beta2*a2)
       qx=qxs+dqn*nxn(mn)     !attention aux signes des normales
       qy=qys+dqn*nyn(mn)
       qz=qzs+dqn*nzn(mn)
!
       v(nl,1)=rho+(pres(m)-ps)*(rho/roc)**2
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
      end subroutine
end module
