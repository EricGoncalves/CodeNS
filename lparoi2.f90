module mod_lparoi2
implicit none
contains
      subroutine lparoi2( &
                 lm,ncyc, &
                 v,mu,mut,dist, &
                 nxn,nyn,nzn, &
                 ncin,ncbd,mfbm, &
                 toxx,toxy,toxz,toyy,toyz,tozz, &
                 qcx,qcy,qcz, &
                 mnpar,fgam,tp, &
                 temp)
!
!***********************************************************************
!
!_DA  DATE_C :fevrier 1999 -- AUTEUR :  Eric Goncalves
!
!     ACT
!_A    Lois de paroi
!_A    - parois isothermes
!_A    - contrainte de frottement imposee a la paroi
!_A    - flux de chaleur imposee a la paroi
!_A
!_A    A partir du champ aerodynamique au niveau des cellules adjacentes
!_A    aux parois, on calcule les densites de flux numeriques a la paroi.
!
!
!     INP
!_I    ncin       : arg int (ip41      ) ; ind dans un tab tous domaines de la
!_I                                        cell. interieure adjacente a la front
!_I    ncbd       : arg int (ip41      ) ; ind dans un tab tous domaines d'une
!_I                                        cellule frontiere fictive
!_I    npfb       : com int (lt        ) ; pointeur fin de domaine precedent
!_I                                        dans tableau toutes facettes
!_I    nnn        : com int (lt        ) ; nombre de noeuds du domaine (dont fictif)
!_I    nxn        : arg real(ip42      ) ; composante en x du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nyn        : arg real(ip42      ) ; composante en y du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    nzn        : arg real(ip42      ) ; composante en z du vecteur directeur
!_I                                        normal a une facette frontiere
!_I    v          : arg real(ip11,ip60 ) ; variables a l'instant n+alpha
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
!_I    cl         : com char(mtb       ) ; type de cond lim a appliquer
!_I    img        : com int              ; niveau de grille (multigrille)
!_I    mtbx       : com int              ; nbr total de frontieres
!_I    mmb        : com int (mtt       ) ; nombre de facettes d'une frontiere
!_I    mpb        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tableaux de base des front.
!_I    nba        : com int (mtb       ) ; rang de traitement d'une front
!_I    ndlb       : com int (mtb       ) ; numero dom contenant la frontiere
!_I    npfb       : com int (lt        ) ; pointeur fin de domaine precedent
!_I                                        dans tableau toutes facettes
!_I    mpn        : com int (mtt       ) ; pointeur fin de front precedente
!_I                                        dans tab front a normales stockees
!_I    nbd        : com int              ; nombre de frontieres a traiter
!_I    lbd        : com int (mtt       ) ; numero de front a traiter
!_I    mnpar      : arg real(ip12      ) ; pointeur dans tableaux front normales
!_I                                        stockees du point de rattach normale
!_I    fgam       : arg real(ip42      ) ; fonction d'intermittence pour
!_I                                        transition
!
!
!-----parameters figes--------------------------------------------------
!
      use para_var
      use para_fige
      use maillage
      use boundary
      use proprieteflu
      use modeleturb
      use definition
implicit none
integer :: lm
integer :: ncyc
double precision :: v
double precision :: dist
integer :: ncin
integer :: ncbd
integer :: mfbm
double precision :: toxx
double precision :: toxy
double precision :: toxz
double precision :: toyy
double precision :: toyz
double precision :: tozz
double precision :: qcx
double precision :: qcy
double precision :: qcz
integer :: mnpar
double precision :: fgam
double precision :: tp
double precision :: temp
double precision :: cta
double precision :: ctb
double precision :: denom
integer :: iter
integer :: m
integer :: m0ns
integer :: mb
integer :: mpar
integer :: mt
integer :: n0c
integer :: nc
integer :: nfacns
integer :: ni
integer :: nii
double precision :: phip
double precision :: qc1
double precision :: rop
double precision :: t1
double precision :: t2
double precision :: t3
double precision :: temp1
double precision :: tn
double precision :: top
double precision :: tt
double precision :: upyp1
double precision :: v1t
double precision :: v1x
double precision :: v1y
double precision :: v1z
double precision :: yp02
!
!-----------------------------------------------------------------------
!
      logical lamin
      real mu,mut,mup,sv
      real nxn,nyn,nzn,n1,n2,n3
!
      dimension mu(ip12),mut(ip12), &
                toxx(ip12),toxy(ip12),toxz(ip12), &
                toyy(ip12),toyz(ip12),tozz(ip12), &
                qcx(ip12),qcy(ip12),qcz(ip12)
      dimension nxn(ip42),nyn(ip42),nzn(ip42)
      dimension ncin(ip41),ncbd(ip41)
      dimension v(ip11,ip60),dist(ip12)
      dimension mnpar(ip12),fgam(ip42)
      dimension tp(ip40)
      dimension temp(ip11)
!
      sv=110.4/tnz !air
      mt=mmb(mfbm)
      m0ns=mpn(mfbm)
      n0c=npc(lm)
!
!     boucle sur les facettes d'une frontiere paroi
      do m=1,mt
        mb=mpb(mfbm)+m
        ni=ncin(mb)
        nc=ncbd(mb)
        nfacns=m0ns+m
        mpar=mnpar(ni)
        nii=ni-n0c
!       test sur transition et regime d'ecoulement
        if((fgam(mpar).lt.1.e-3).and.(ktransi.gt.0)) then
!         laminaire
          lamin=.true.
        else
!         turbulent
          lamin=.false.
        endif
!       vitesse cellule adjacente a la paroi (cellule 1)
        v1x=v(ni,2)/v(ni,1)
        v1y=v(ni,3)/v(ni,1)
        v1z=v(ni,4)/v(ni,1)
!       normale a la paroi
        n1=nxn(nfacns)
        n2=nyn(nfacns)
        n3=nzn(nfacns)
!       tangente normee a la paroi
        tn=v1x*n1+v1y*n2+v1z*n3
        t1=v1x-tn*n1
        t2=v1y-tn*n2
        t3=v1z-tn*n3
        tt=sqrt(t1**2+t2**2+t3**2)
        t1=t1/tt
        t2=t2/tt
        t3=t3/tt
!       composante tangentielle de la vitesse dans repere paroi : v1t
        v1t=v1x*t1+v1y*t2+v1z*t3
!       temperature cellule 1 : temp1
        temp1=temp(ni)
!       masse volumique a la paroi
        rop=v(ni,1)*temp1/tp(m)
!       viscosite moleculaire a la paroi
        mup=mu(ni)*sqrt(tp(m)/temp1)*(1.+sv/temp1)/(1.+sv/tp(m))
!       correction de compressibilite (loi de Van Driest)
        cta=(mu(ni)+mut(ni))/(cp*(mu(ni)/pr+mut(ni)/prt))
        ctb=cta/(2.*tp(m))
        denom=(tp(m)-temp1-cta*0.5*v1t**2)*sqrt(temp1/tp(m)) + &
               tp(m)-temp1+cta*0.5*v1t**2
        v1t=(1./sqrt(ctb))*asin(2.*sqrt(ctb)*(tp(m)-temp1)*v1t/denom)
!       contrainte de frottement a la paroi : top
        upyp1=rop*v1t*dist(ni)/mup
        yp02=yp0**2
        if(upyp1.le.yp02 .or. lamin) then
!         loi lineaire
          top=mup*v1t/dist(ni)
        else
!         loi logarithmique
          top=mup*v1t/dist(ni)
          do iter=1,10
            top=rop*v1t**2/(log(dist(ni)*sqrt(rop*top)/mup)/vkar+cllog)**2
          enddo
        endif
!
!         if(lamin) then
!         loi lineaire
!          top=mup*v1t/dist(ni)
!         else
!         loi de Spalding
!          top=mup*v1t/dist(ni)
!          upl=v1t*sqrt(rop/top)
!          do jj=1,15
!            fu=rop*dist(ni)*v1t/(mup*upl) - upl - exp(-vkar*cllog)*
!     &         (exp(vkar*upl) -1. -vkar*upl - 0.5*(vkar*upl)**2 -
!     &          ((vkar*upl)**3)/6.)
!            dfu=-rop*dist(ni)*v1t/(mup*upl**2)-1.-exp(-vkar*cllog)*
!     &          vkar*(exp(vkar*upl) -1. -vkar*upl - 0.5*(vkar*upl)**2)
!            upl=upl-fu/dfu
!          enddo
!         endif
!          top=rop*(v1t/upl)**2
!
!        loi de paroi raffine (interpolation dans la zone tampom)
!        if(upyp1.le.9.) then
!          loi lineaire
!           top=mup*v1t/dist(ni)
!       else
!         yplus0=3.
!         do ii=1,10
!           yplus0=upyp1/(log(yplus0)/vkar+cllog)
!          end do
!           if(yplus0.gt.40.) then
!            loi logarithmique
!           top=mup*v1t/dist(ni)
!           do iter=1,10
!            top=rop*v1t**2/(log(dist(ni)*sqrt(rop*top)/mup)/
!     &               vkar+cllog)**2
!            end do
!           else
!           region tampon : interpolation avec polynome de degre 4
!           top=mup*v1t/dist(ni)
!           do jj=1,10
!             yy=log(dist(ni)*sqrt(rop*top)/mup)
!               top=max(1.e-10,rop*v1t** 2/(0.17962*yy**4-
!     &                  2.2117*yy**3+9.2052*yy**2-10.804*yy
!     &             +6.4424)**2)
!            end do
!          end if
!      end if
!
!       flux de chaleur a la paroi
        phip=top*(tp(m)-temp1-0.5*cta*v1t**2)/(cta*v1t)
!       tenseur des contraintes a la paroi dans repere general
        toxx(nc)=2*t1*n1*top
        toyy(nc)=2*t2*n2*top
        tozz(nc)=2*t3*n3*top
        toxy(nc)=(t2*n1+t1*n2)*top
        toxz(nc)=(t1*n3+t3*n1)*top
        toyz(nc)=(t2*n3+t3*n2)*top
!       tenseur des contraintes en 1 dans repere general
        toxx(ni)=toxx(nc)
        toyy(ni)=toyy(nc)
        tozz(ni)=tozz(nc)
        toxy(ni)=toxy(nc)
        toxz(ni)=toxz(nc)
        toyz(ni)=toyz(nc)
!       flux de chaleur a la paroi dans repere general
!       ATTENTION! le signe moins provient de la convention utilisee dans le code
        qcx(nc)=-n1*phip
        qcy(nc)=-n2*phip
        qcz(nc)=-n3*phip
!       flux de chaleur dans cellule 1 dans repere general
!       ATTENTION! le signe moins provient de la convention utilisee dans le code
        qc1=v1t*top+phip
        qcx(ni)=-n1*qc1
        qcy(ni)=-n2*qc1
        qcz(ni)=-n3*qc1
!
!     fin boucle sur facettes d'une frontiere paroi
      enddo
!
      return
      end
end module
