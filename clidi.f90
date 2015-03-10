      subroutine clidi( &
                 mfb,l,rpi,rti, &
                 nxn,nyn,nzn,ncbd,v, &
                 mmb,mpb,mpn, &
                 pression,temp,cson)
!
!***********************************************************************
!
!_DA  DATE_C : decembre 2007  - Eric GONCALVES / LEGI
!
!     ACT
!_A    Calcul des variables sur les facettes frontieres par
!_A    traitement des conditions aux limites :
!_A     pression d'arret imposee , pi=pa1*rpi
!_A     enthalpie d'arret imposee ,  hi=ha1*rti
!       normales interieures.
!
!     INP
!_I    mfb        : arg int              ; numero de frontiere
!_I    l          : arg int              ; numero de domaine
!_I    rpi        : arg real(ip40      ) ; pres d'arret/pres d'arret etat de
!_I                                        ref utilisateur a imposer
!_I    rti        : arg real(ip40      ) ; temp d'arret/temp d'arret etat de
!_I                                        ref utilisateur a imposer
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n
!_I    mmb        : arg int (mtb       ) ; nombre de pts d'une frontiere
!_I    mpb        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    mpn        : arg int (mtb       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    gam        : com real             ; rapport des chaleurs specifiques
!_I    gam1       : com real             ; rap chal spec -1
!_I    gam2       : com real             ; (rap chal spec -1)/2
!_I    gam4       : com real             ; 1/(rap chal spec -1)
!_I    gam5       : com real             ; rap chal spec/(rap chal spec -1)
!_I    aa1        : com real             ; etat de reference utilisateur
!_I                                        adimensionne, vitesse du son d'arret
!_I    pa1        : com real             ; pression d'arret de l'etat
!_I                                        de reference utilisateur adimensionne
!_I    id1        : com int (lt        ) ; indice min en i fictif
!_I    id2        : com int (lt        ) ; indice max en i fictif
!_I    jd1        : com int (lt        ) ; indice min en j fictif
!_I    jd2        : com int (lt        ) ; indice max en j fictif
!_I    npn        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab tous noeuds
!_I    npc        : com int (lt        ) ; pointeur fin de dom precedent
!_I                                        dans tab toutes cellules
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
      use schemanum
      use definition
implicit none
integer :: inc
integer :: mfb
integer :: l
double precision :: rpi
double precision :: rti
integer :: ncbd
double precision :: v
integer :: mmb
integer :: mpb
integer :: mpn
double precision :: pression
double precision :: temp
double precision :: cson
integer :: id
integer :: jd
integer :: kd
double precision :: b
double precision :: gam2t
double precision :: gam6
double precision :: gamt
integer :: m
integer :: mb
integer :: mn
integer :: mt
integer :: n0c
integer :: n0n
integer :: nc
integer :: nci
integer :: ncij
integer :: ncijk
integer :: ncik
integer :: ncj
integer :: ncjk
integer :: nck
integer :: nid
integer :: nijd
integer :: njd
double precision :: pres
double precision :: qinf
double precision :: qn
double precision :: qxs
double precision :: qys
double precision :: qzs
double precision :: rho
!
!-----------------------------------------------------------------------
!
      real nxn,nyn,nzn
!
      dimension rpi(ip40),rti(ip40)
      dimension v(ip11,ip60),pression(ip11),temp(ip11),cson(ip11)
      dimension nxn(ip42),nyn(ip42),nzn(ip42),ncbd(ip41)
      dimension mmb(mtt),mpb(mtt),mpn(mtt)
!
      inc(id,jd,kd)=id+jd*nid+kd*nijd
!
      n0n=npn(l)
      n0c=npc(l)
!
      nid = id2(l)-id1(l)+1
      njd = jd2(l)-jd1(l)+1
      nijd = nid*njd
!
      nci = inc(1,0,0)
      ncj = inc(0,1,0)
      nck = inc(0,0,1)
      ncij = inc(1,1,0)
      ncik = inc(1,0,1)
      ncjk = inc(0,1,1)
      ncijk= inc(1,1,1)
!
      mt=mmb(mfb)
      gam6=gam/gam1
      gamt=gam/aa1**2
      gam2t=0.5*gam1/aa1**2
      qinf=rm0*aa1/(1.+gam2*rm0**2)**0.5
!
      do m=1,mt
       mb=mpb(mfb)+m
       mn=mpn(mfb)+m
       nc=ncbd(mb)
!
       qxs=v(nc,2)/v(nc,1)
       qys=v(nc,3)/v(nc,1)
       qzs=v(nc,4)/v(nc,1)
       qn=qxs*nxn(mn)+qys*nyn(mn)+qzs*nzn(mn)
!-----la notation b designe la quantite 1/(1+((gam-1)/2)*mach**2)
!       b=1./(1.+gam2t/rti(m)*qn**2)
       b=1./(1.+gam2t/rti(m)*qinf**2)
       pres=pa1*rpi(m)*b**gam6
       rho=gamt*pres/(b*rti(m))
!
       v(nc,1)=rho
!       v(nc,2)=rho*qn
       v(nc,2)=rho*qinf
       v(nc,3)=0.
       v(nc,4)=0.
       v(nc,5)=pres/gam1+pinfl+0.5*(v(nc,2)**2)/rho
!
       pression(nc)=pres
       temp(nc)=gam*pres/rho
       cson(nc)=sqrt(temp(nc))
      enddo
!
      return
      end
