module mod_cldebit_prcd
implicit none
contains
      subroutine cldebit_prcd( &
                 mfb,pres, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn,l, &
                 pression,temp,cson)
!
!***********************************************************************
!
!_DA           janvier 2006 - Eric GONCALVES / LEGI
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement de la condition aux limites : debit impose.
!_A    Utilisation des quatre relations de compatibilite discretisees
!_A    associee aux quatre caracteristiques sortantes
!_A      Dut=0
!_A      Dut=0
!_A      Dp - c2*Drho = 0
!_A      (lambda+ -un)*Dp + rho*beta2*c2*Dun = 0    (vitesse normale)
!_A    Normales interieures.
!_A    ATTENTION: le debit local est sotcke dans le tableau 'pres'.
!_A    Precondtionnement basse vitesse de Turkel
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
      use schemanum 
      use definition
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
double precision :: q2
double precision :: qinf
double precision :: qns
double precision :: qx
double precision :: qxs
double precision :: qy
double precision :: qys
double precision :: qz
double precision :: qzs
double precision :: rhos
!
!-----------------------------------------------------------------------
!
      real nxn,nyn,nzn
      dimension pres(ip40)
      dimension v(ip11,ip60)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
      dimension mmb(mtt),mpb(mtt),mpn(mtt)
      dimension pression(ip11),temp(ip11),cson(ip11)
!
      qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
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
       q2=qxs**2+qys**2+qzs**2
       ps=gam1*(v(nl,5)-0.5*rhos*q2-pinfl)
       a2=gam*ps/rhos
       beta2=min(max(q2/a2,cte*qinf**2/a2),1.)
       alp=0.5*((1.+beta2)*abs(qns) + sqrt(((1.-beta2)*qns)**2+4.*beta2*a2))
!
       ca=(alp+qns)*a2
       cb=-ca*rhos+rhos*beta2*a2*qns
       cc=rhos*beta2*a2*pres(m)
       drho=0.5*(-cb+sqrt(abs(cb**2-4.*ca*cc)))/ca-rhos
       dp=a2*drho
       dqn=(alp+qns)*dp/(rhos*beta2*a2)
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
